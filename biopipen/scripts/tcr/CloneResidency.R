{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}

library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(glue)
library(ggnewscale)
library(ggplot2)
library(ggprism)
library(ggVennDiagram)
library(ComplexUpset)

theme_set(theme_prism())


immfile <- {{ in.immdata | r }}
metafile <- {{ in.metafile | r }}
outdir <- {{ out.outdir | r }}
joboutdir <- {{ job.outdir | r }}

subject_key <- {{ envs.subject | r }}
group_key <- {{ envs.group | r }}
sample_order <- {{ envs.order | r }}
section <- {{ envs.section | r }}
mutaters <- {{ envs.mutaters | r }}
subset <- {{ envs.subset | r }}
prefix <- {{ envs.prefix | r }}
upset_ymax <- {{ envs.upset_ymax | r }}
upset_trans <- {{ envs.upset_trans | r }}
cases <- {{ envs.cases | r }}

# Fill up cases using `envs.xxx` if not provided and compose a DEFAULT case
# if no cases are provided
log_info("Preparing cases...")
if (is.null(cases) || length(cases) == 0) {
    cases <- list(
        DEFAULT = list(
            subject = subject_key,
            group = group_key,
            order = sample_order,
            subset = subset,
            section = section,
            upset_ymax = upset_ymax,
            upset_trans = upset_trans
        )
    )
} else {
    for (key in names(cases)) {
        cases[[key]]$subject <- cases[[key]]$subject %||% subject_key
        cases[[key]]$group <- cases[[key]]$group %||% group_key
        cases[[key]]$order <- cases[[key]]$order %||% sample_order
        cases[[key]]$section <- cases[[key]]$section %||% section
        cases[[key]]$subset <- cases[[key]]$subset %||% subset
        cases[[key]]$upset_ymax <- cases[[key]]$upset_ymax %||% upset_ymax
        cases[[key]]$upset_trans <- cases[[key]]$upset_trans %||% upset_trans
    }
}

log_info("Preparing data ...")
log_info("- Loading immdata ...")
immdata <- readRDS(immfile)

log_info("- Expanding immdata$data to cell-level ...")
cldata <- do_call(rbind, lapply(names(immdata$data), function(name) {
    dat <- immdata$data[[name]] %>% separate_rows(Barcode, sep = ";") # Split barcodes
    dat$Sample <- name
    dat <- dat %>% left_join(immdata$meta, by = "Sample", suffix = c("", "_meta"))

    if (!is.null(prefix) && nchar(prefix) > 0) {
        dat <- dat %>% mutate(Barcode = glue(paste0(prefix, "{Barcode}")))
    }
}))

if (!is.null(metafile)) {
    log_info("- Loading metafile ...")
    # Check if extension is rds/RDS, if so, it should be a Seurat object
    if (endsWith(metafile, ".rds") || endsWith(metafile, ".RDS")) {
        meta <- readRDS(metafile)@meta.data
    } else {
        meta <- read.table(
            metafile, row.names = 1, sep = "\t", header = TRUE, stringsAsFactors = FALSE
        )
    }

    log_info("- Merging metafile to cldata ...")
    cldata <- cbind(
        cldata,
        meta[cldata$Barcode, setdiff(colnames(meta), colnames(cldata)), drop = FALSE]
    )
}

log_info("Applying mutaters ...")
if (!is.null(mutaters) && length(mutaters) > 0) {
    cldata <- cldata %>% mutate(!!!lapply(mutaters, parse_expr))
}

# Scatter plot functions
exponent <- function(x) {
    floor(log10(abs(x)))
}

mantissa <- function(x) {
    mant <- log10(abs(x))
    10 ^ (mant - floor(mant))
}

get_groups <- function(order) {
    # order is something like [`A,B`, `A,C`, `B,C`]
    # return `A`, `B`, `C`
    unique(unlist(strsplit(order, ",")))
}

perpare_case <- function(casename, case) {
    log_info("- Processing case: {casename} ...")
    # Check if required keys are provided
    if (is.null(case$subject) || length(case$subject) == 0) {
        stop(paste("  `subject` is required for case:", casename))
    }
    if (is.null(case$group) || length(case$group) == 0) {
        stop(paste("  `group` is required for case:", casename))
    }
    if (!is.null(case$order) && length(case$order) > 0) {
        has_comma <- grepl(",", case$order)
        if (all(has_comma)) {
            # It's recommended
            case$order <- unname(sapply(
                case$order,
                function(x) paste(trimws(strsplit(x, ",")[[1]]), collapse=",")
            ))
        } else if (!any(has_comma)) {
            if (length(case$order) > 2) {
                log_warn("  Order of groups is not recommended, please use comma to separate groups.")
                log_warn("  Instead of `['A', 'B', 'C']`, use `['A,B', 'A,C', 'B,C']`.")
                case$order <- sapply(
                    combn(case$order, 2, simplify = FALSE),
                    function(x) paste(x, collapse = ",")
                )
            } else {
                case$order <- paste(case$order, collapse = ",")
            }
        } else {
            stop(
                paste0(
                    "  Order of groups in case:", casename,
                    " is not consistent, please use comma to separate groups. ",
                    "Instead of `['A', 'B', 'C']`, use `['A,B', 'A,C', 'B,C']`, ",
                    "however, this is inconsistent: `['A,B', 'C']`"
                )
            )
        }
    }

    # Create case-specific directories
    # Scatter plots
    scatter_dir <- file.path(outdir, slugify(casename), "scatter")
    dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)

    # Venn diagrams
    venn_dir <- file.path(outdir, slugify(casename), "venn")
    dir.create(venn_dir, recursive = TRUE, showWarnings = FALSE)

    # Upset plots
    upset_dir <- file.path(outdir, slugify(casename), "upset")
    dir.create(upset_dir, recursive = TRUE, showWarnings = FALSE)

    # Counts
    counts_dir <- file.path(outdir, slugify(casename), "counts")
    dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)

    case
}

plot_scatter <- function(counts, subject, suf1, suf2) {
    # options(repr.plot.width=9, repr.plot.height=7)

    dual <- which(counts[[suf1]] > 0 & counts[[suf2]] > 0)
    if (length(dual) <= 2) {
        test <- list(estimate = NA, p.value = NA)
    } else {
        test <- cor.test(log(counts[[suf1]][dual]), log(counts[[suf2]][dual]))
    }
    sum_counts1 <- sum(counts[[suf1]])
    sum_counts2 <- sum(counts[[suf2]])

    counts1_norm <- jitter(1 + counts[[suf1]], amount = 0.25) / sum_counts1
    counts2_norm <- jitter(1 + counts[[suf2]], amount = 0.25) / sum_counts2

    oo <- sample(length(counts1_norm))
    plotdata <- data.frame(x = counts1_norm[oo], y = counts2_norm[oo])
    # plotdata$color = cl.colors[oo]
    names(plotdata) <- c(suf1, suf2)
    plotdata <- plotdata %>% mutate(
        Type = case_when(
            counts[[suf1]][oo] == 1 & counts[[suf2]][oo] == 0 ~ paste(suf1, "Singleton"),
            counts[[suf1]][oo] == 0 & counts[[suf2]][oo] == 1 ~ paste(suf2, "Singleton"),
            counts[[suf1]][oo] > 1 & counts[[suf2]][oo] == 0 ~ paste(suf1, "Multiplet"),
            counts[[suf1]][oo] == 0 & counts[[suf2]][oo] > 1 ~ paste(suf2, "Multiplet"),
            counts[[suf1]][oo] < counts[[suf2]][oo] ~ "Expanded",
            counts[[suf1]][oo] > counts[[suf2]][oo] ~ "Collapsed",
            TRUE ~ "Dual"
        ),
        Type = as.factor(Type),
        Size = pmax(counts1_norm[oo], counts2_norm[oo])
    )

    xbreaks <- c(
        1 / sum_counts1,
        0.001 + 1 / sum_counts1,
        0.01 + 1 / sum_counts1,
        0.1 + 1 / sum_counts1
    )
    ybreaks <- c(
        1 / sum_counts2,
        0.001 + 1 / sum_counts2,
        0.01 + 1 / sum_counts2,
        0.1 + 1 / sum_counts2
    )

    minx <- min(plotdata[[suf1]])
    miny <- min(plotdata[[suf2]])
    maxx <- max(plotdata[[suf1]])
    maxy <- max(plotdata[[suf2]])
    # color = plotdata$color
    # names(color) = color
    # patient = as.character(patient)
    n_formatted <- formatC(length(oo), format = "f", big.mark = ",", digits = 0)
    r_formatted <- format(test$estimate, digits = 2, scientific = F)
    if (is.na(test$p.value)) {
        subtitle <- bquote(
            italic(n)[D] == .(length(dual)) ~ ~ italic(r) == .(r_formatted) ~ ~ italic(P) == "NA"
        )
    } else if (test$p.value < 1e-4) {
        P_mant <- format(mantissa(test$p.value), digits = 2)
        P_exp <- exponent(test$p.value)
        subtitle <- bquote(
            italic(n)[D] == .(length(dual)) ~ ~ italic(r) ==
            .(r_formatted) ~ ~ italic(P) == .(P_mant) %*% 10^.(P_exp)
        )
    } else {
        P_formatted <- format(test$p.value, digits = 2)
        subtitle <- bquote(
            italic(n)[D] == .(length(dual)) ~ ~ italic(r) ==
            .(r_formatted) ~ ~ italic(P) == .(P_formatted)
        )
    }
    ggplot(plotdata) +
        geom_point(
            aes(
                x = !!sym(suf1),
                y = !!sym(suf2),
                color = Type,
                size = Size,
                fill = Type
            ),
            alpha = .6,
            shape = 21
        ) +
        scale_x_continuous(
            trans = "log2",
            limits = c(minx, maxx),
            breaks = xbreaks,
            labels = c("0", "0.001", "0.01", "0.1")
        ) +
        scale_y_continuous(
            trans = "log2",
            limits = c(miny, maxy),
            breaks = ybreaks,
            labels = c("0", "0.001", "0.01", "0.1")
        ) +
        theme_prism(base_size = 16) +
        scale_size(guide = "none") +
        labs(
            title = bquote(.(subject) ~ (italic(n) == .(n_formatted))),
            subtitle = subtitle
        ) +
        geom_segment(
            data = data.frame(
                # diagnal, horizontal, vertical, horizontal short, vertical short
                x = c(1.5 / sum_counts1, minx, 1.5 / sum_counts1, minx, 2.5 / sum_counts1),
                xend = c(maxx, maxx, 1.5 / sum_counts1, 1.5 / sum_counts1, 2.5 / sum_counts1),
                y = c(1.5 / sum_counts2, 1.5 / sum_counts2, miny, 2.5 / sum_counts2, miny),
                yend = c(maxy, 1.5 / sum_counts2, maxy, 2.5 / sum_counts2, 1.5 / sum_counts2)
            ),
            aes(x = x, y = y, xend = xend, yend = yend), color = "gray"
        )
}

plot_venndg <- function(counts, groups, singletons) {
    venn_data <- list()
    for (group in groups) {
        venn_data[[group]] <- counts %>% filter(!!sym(group) > 0) %>% pull(CDR3.aa)
    }
    venn <- Venn(venn_data)
    vdata <- process_data(venn)
    vregion <- venn_region(vdata)
    vregion$singleton_count = singletons[vregion$name, "count"]
    vregion <- vregion %>% mutate(
        count_perc = round(count / sum(count) * 100, 1),
        count_str = paste0(count, " (", count_perc, "%)"),
        count_str = if_else(is.na(singleton_count), count_str, paste0(count_str, "\nsingletons = ", singleton_count))
    )

    venn_p <- ggplot() +
        # 1. region count layer
        geom_sf(aes(fill = count), data = venn_region(vdata)) +
        # 2. set edge layer
        # geom_sf(aes(color = factor(id)), data = venn_setedge(data), show.legend = FALSE) +
        # 3. set label layer
        geom_sf_text(aes(label = name), data = venn_setlabel(vdata)) +
        # 4. region label layer
        geom_sf_label(
            aes(label = count_str),
            alpha = .8,
            label.padding = unit(.2, "lines"),
            data = vregion
        ) +
        # 5. singletons label layer
        scale_fill_distiller(palette = "Oranges", direction = 1) +
        theme_void() +
        theme(plot.margin = margin(1,1,1,1, "cm"))

    venn_p
}

plot_upset <- function(counts, singletons, upset_ymax, upset_trans) {

    cnts <- column_to_rownames(counts, "CDR3.aa") %>%
        mutate(across(everything(), ~ as.integer(as.logical(.x))))
    sgltns <- unlist(singletons$CDR3.aa)
    cnts$..type <- "Multiplets"
    cnts[sgltns, "..type"] <- "Singletons"
    sets <- setdiff(colnames(cnts), "..type")

    p <- ggplot(mapping = aes(x = intersection, fill = ..type)) +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            panel.grid = element_blank(),
            axis.line.y = element_line(color = "#3b3b3b"),
            axis.ticks.y = element_line(color = "#3b3b3b"),
        ) +
        scale_fill_manual(values = c("#3b3b3b", "orange")) +
        ylab("Intersection size")

    if (is.null(upset_trans)) {
        p <- p + geom_bar(stat = "count", position = "stack") +
            geom_text(
                aes(label = ..count.., vjust = ifelse(..type == "Multiplets", -0.25, +1.25)),
                stat = "count", position = "stack", size = 2.8)
        if (!is.null(upset_ymax)) {
            p <- p + ylim(0, upset_ymax)
        }
    } else {
        p <- p + geom_bar(stat = "count", position = "dodge2") +
            geom_text(
                aes(label = ..count..),
                stat = "count", position = position_dodge(width = 0.9), vjust = -0.25, size = 2.5)

        # limit the y and do log10 transformation
        if (!is.null(upset_ymax)) {
            p <- p + scale_y_continuous(trans = "log10", limits = c(1, upset_ymax))
        } else {
            p <- p + scale_y_continuous(trans = "log10")
        }
    }

    upset(
        cnts, rev(sets),
        sort_sets = FALSE,
        # Remove the base annotations
        base_annotations = list(),
        annotations = list('IntersectSize' = p),
        themes = upset_modify_themes(
            list(
                intersections_matrix = theme(
                    axis.line = element_line(color = "#3b3b3b"),
                    axis.ticks.y = element_line(color = "#3b3b3b"),
                )
            )
        )
    )
}

headings <- function(section, casename, subject) {
    list(h1 = ifelse(
            is.null(section),
            ifelse(casename == "DEFAULT", subject, casename),
            section
        ),
        h2 = ifelse(
            is.null(section),
            ifelse(casename == "DEFAULT", "#", subject),
            ifelse(casename == "DEFAULT", subject, casename)
        ),
        h3 = ifelse(
            is.null(section),
            "#",
            ifelse(casename == "DEFAULT", "#", subject)
        )
    )
}

handle_subject <- function(i, subjects, casename, case) {
    casedir = file.path(outdir, slugify(casename))
    # Generate a residency table
    # |    CDR3.aa    | Tumor | Normal |
    # | SEABESRWEFAEF | 0     | 10     |
    # | AWEARWGAWGGGR | 21    | 1      |
    # | GREWFQQWFEWF  | 34    | 0      |
    subject_row <- subjects[i, , drop = FALSE]
    subject <- subject_row %>%
        select(all_of(case$subject)) %>%
        mutate(across(everything(), as.character)) %>%
        paste(collapse = "-")

    log_info("  Handling {i}/{nrow(subjects)}: {subject} ...")

    if (!is.null(case$subset)) {
        counts <- cldata %>% filter(!!parse_expr(case$subset))
    } else {
        counts <- cldata
    }
    # CDR3.aa | Group | .n
    # --------+-------+----
    # AAAAAAA | Tumor | 10
    # AAAAAAA | Normal| 20
    counts <- subject_row %>%
        left_join(counts, by = case$subject) %>%
        group_by(CDR3.aa, !!sym(case$group)) %>%
        summarise(.n = n())

    if (!is.null(case$order)) {
        groups <- get_groups(case$order)
        counts <- counts %>% filter(!!sym(case$group) %in% groups)
        groups <- intersect(groups, unique(counts[[case$group]]))
    } else {
        groups <- counts %>% pull(!!sym(case$group)) %>% unique()
        case$order <- sapply(combn(groups, 2, simplify = FALSE), function(x) paste(x, collapse = ","))
    }
    if (length(unique(counts[[case$group]])) < 2) {
        log_warn("  - Subject doesn't have enough groups: {subject}")
        return()
    }
    singletons = counts %>%
        group_by(CDR3.aa) %>%
        summarise(name = if_else(sum(.n) == 1, paste(!!sym(case$group), collapse=""), "")) %>%
        filter(name != "") %>%
        group_by(name) %>%
        summarise(count = length(CDR3.aa), CDR3.aa = list(CDR3.aa)) %>%
        column_to_rownames("name")

    counts <- counts %>%
        pivot_wider(
            id_cols = CDR3.aa,
            names_from = !!sym(case$group),
            values_from = .n
        ) %>%
        select(CDR3.aa, !!!syms(groups))
    counts[is.na(counts)] <- 0

    # Save counts
    counts_dir <- file.path(casedir, "counts")
    countfile <- file.path(counts_dir, paste0(slugify(subject), ".txt"))
    write.table(
        counts,
        file = countfile,
        sep = "\t",
        row.names = TRUE,
        col.names = TRUE,
        quote = FALSE
    )

    h <- headings(case$section, casename, "Clone Size Tables")
    add_report(
        list(kind = "table", src = countfile, ds_name = subject),
        h1 = h$h1,
        h2 = h$h2,
        h3 = h$h3,
        ui = "dropdown_switcher"
    )

    # scatter plot
    # Make plots B ~ A, C ~ B, and C ~ A for order A, B, C
    # combns <- combn(groups, 2, simplify = FALSE)
    h <- headings(case$section, casename, "Residency Plots")
    scatter_dir <- file.path(casedir, "scatter")
    for (j in seq_along(case$order)) {
        pair <- strsplit(case$order[j], ",")[[1]]
        if (length(setdiff(pair, groups)) > 0) {
            log_warn("  - Comparison {case$order[j]} doesn't exist.")
            next
        }
        scatter_p <- plot_scatter(counts, subject, pair[1], pair[2])
        scatter_png <- file.path(
            scatter_dir,
            paste0("scatter_", slugify(subject), "_", slugify(pair[1]), "_", slugify(pair[2]), ".png")
        )
        png(scatter_png, res = 100, height = 800, width = 1000)
        print(scatter_p)
        dev.off()

        scatter_pdf <- gsub(".png$", ".pdf", scatter_png)
        pdf(scatter_pdf, width = 10, height = 8)
        print(scatter_p)
        dev.off()

        add_report(
            list(
                name = paste0(subject, " (", pair[1], " - ", pair[2], ")"),
                src = scatter_png,
                download = scatter_pdf
            ),
            h1 = h$h1,
            h2 = h$h2,
            h3 = h$h3,
            ui = "table_of_images:3"
        )
    }

    # upset/venn
    venn_dir <- file.path(casedir, "venn")
    venn_png <- file.path(venn_dir, paste0("venn_", slugify(subject), ".png"))
    venn_pdf <- gsub(".png$", ".pdf", venn_png)
    p <- plot_venndg(counts, groups, singletons)
    png(venn_png, res = 100, height = 600, width = 800)
    print(p)
    dev.off()

    pdf(venn_pdf, width = 8, height = 6)
    print(p)
    dev.off()

    h <- headings(case$section, casename, "Overlapping Clones (Venn Diagram)")
    add_report(
        list(src = venn_png, name = subject, download = venn_pdf),
        h1 = h$h1,
        h2 = h$h2,
        h3 = h$h3,
        ui = "table_of_images:3"
    )

    upset_dir <- file.path(casedir, "upset")
    upset_png <- file.path(upset_dir, paste0("upset_", slugify(subject), ".png"))
    upset_pdf <- gsub(".png$", ".pdf", upset_png)
    p <- plot_upset(counts, singletons, case$upset_ymax, case$upset_trans)
    png(upset_png, res = 100, height = 600, width = 800)
    print(p)
    dev.off()

    pdf(upset_pdf, width = 8, height = 6)
    print(p)
    dev.off()

    h <- headings(case$section, casename, "Overlapping Clones (UpSet Plots)")
    add_report(
        list(src = upset_png, name = subject, download = upset_pdf),
        h1 = h$h1,
        h2 = h$h2,
        h3 = h$h3,
        ui = "table_of_images:3"
    )
}

handle_case <- function(casename, case) {
    case <- perpare_case(casename, case)
    if (!is.null(case$subset)) {
        subjects <- cldata %>%
            filter(!!parse_expr(case$subset)) %>%
            distinct(!!!syms(case$subject)) %>%
            drop_na()
    } else {
        subjects <- cldata %>%
            distinct(!!!syms(case$subject)) %>%
            drop_na()
    }

    h <- headings(case$section, casename, "Clone Size Tables")
    add_report(
        list(ds_name = "Select a subject ..."),
        h1 = h$h1,
        h2 = h$h2,
        h3 = h$h3,
        ui = "dropdown_switcher"
    )

    h <- headings(case$section, casename, "Residency Plots")
    add_report(
        list(
            kind = "descr",
            content = paste0(
                "The residency plots showing the clones of paired samples (x-axis and y-axis). ",
                "The size of the dot represents the relative abundance of the clone. ",
                "The color of the dot represents the type of the clone: "
            )
        ),
        list(
            kind = "list",
            items = c(
                "Collapsed (clones that are less abundant in the y-axis sample)",
                "Dual (clones that are equally abundant in both samples)",
                "Expanded (clones that are more abundant in the y-axis sample)",
                "(x-axis sample) Multiplet (clones that are only present in the x-axis sample, with multiple cells)",
                "(x-axis sample) Singleton (clones that are only present in the x-axis sample, with a single cell)",
                "(y-axis sample) Multiplet (clones that are only present in the y-axis sample, with multiple cells)",
                "(y-axis sample) Singleton (clones that are only present in the y-axis sample, with a single cell)"
            )
        ),
        h1 = h$h1,
        h2 = h$h2,
        h3 = h$h3,
        ui = "flat"
    )

    h <- headings(case$section, casename, "Overlapping Clones (Venn Diagram)")
    add_report(
        list(
            kind = "descr",
            content = "For samples in each subject, showing the overlapping clones between samples in Venn diagrams."
        ),
        h1 = h$h1,
        h2 = h$h2,
        h3 = h$h3,
        ui = "flat"
    )

    h <- headings(case$section, casename, "Overlapping Clones (UpSet Plots)")
    add_report(
        list(
            kind = "descr",
            content = "For samples in each subject, showing the overlapping clones between samples in UpSet plots."
        ),
        h1 = h$h1,
        h2 = h$h2,
        h3 = h$h3,
        ui = "flat"
    )
    sapply(seq_len(nrow(subjects)), handle_subject, subjects, casename, case)
}

for (casename in sort(names(cases))) {
    handle_case(casename, cases[[casename]])
}

save_report(joboutdir)
