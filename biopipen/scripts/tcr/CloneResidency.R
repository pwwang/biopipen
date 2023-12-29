source("{{biopipen_dir}}/utils/misc.R")

library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(glue)
library(ggnewscale)
library(ggplot2)
library(ggprism)
library(ggVennDiagram)
library(UpSetR)
library(slugify)

theme_set(theme_prism())


immfile <- {{ in.immdata | quote }}
metafile <- {{ in.metafile | r }}
outdir <- {{ out.outdir | quote }}
joboutdir <- {{ job.outdir | quote }}

subject_key <- {{ envs.subject | r }}
group_key <- {{ envs.group | r }}
sample_order <- {{ envs.order | r }}
section <- {{ envs.section | r }}
mutaters <- {{ envs.mutaters | r }}
subset <- {{ envs.subset | r }}
prefix <- {{ envs.prefix | r }}
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
            section = section
        )
    )
} else {
    for (key in names(cases)) {
        if (is.null(cases[[key]]$subject)) {
            cases[[key]]$subject <- subject_key
        }
        if (is.null(cases[[key]]$group)) {
            cases[[key]]$group <- group_key
        }
        if (is.null(cases[[key]]$order)) {
            cases[[key]]$order <- sample_order
        }
        if (is.null(cases[[key]]$section)) {
            cases[[key]]$section <- section
        }
        if (is.null(cases[[key]]$subset)) {
            cases[[key]]$subset <- subset
        }
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
    log_info("Preparing case: {casename} ...")
    # Check if required keys are provided
    if (is.null(case$subject) || length(case$subject) == 0) {
        stop(paste("`subject` is required for case:", casename))
    }
    if (is.null(case$group) || length(case$group) == 0) {
        stop(paste("`group` is required for case:", casename))
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
                log_warn(
                    paste0(
                        "- Order of groups in case:", casename,
                        " is not recommended, please use comma to separate groups. \n",
                        "Instead of `['A', 'B', 'C']`, use `['A,B', 'A,C', 'B,C']`."
                    )
                )
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
                    "- Order of groups in case:", casename,
                    " is not consistent, please use comma to separate groups. \n",
                    "Instead of `['A', 'B', 'C']`, use `['A,B', 'A,C', 'B,C']`, ",
                    "however, this is inconsistent: `['A,B', 'C']`"
                )
            )
        }
    }

    # Create case-specific directories
    # Scatter plots
    scatter_dir <- file.path(outdir, slugify(casename, tolower = FALSE), "scatter")
    dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)

    # Venn diagrams
    venn_dir <- file.path(outdir, slugify(casename, tolower = FALSE), "venn")
    dir.create(venn_dir, recursive = TRUE, showWarnings = FALSE)

    # Upset plots
    upset_dir <- file.path(outdir, slugify(casename, tolower = FALSE), "upset")
    dir.create(upset_dir, recursive = TRUE, showWarnings = FALSE)

    # Counts
    counts_dir <- file.path(outdir, slugify(casename, tolower = FALSE), "counts")
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
            aes_string(
                x = bQuote(suf1), y = bQuote(suf2), color = "Type", size = "Size", fill = "Type"
            ),
            alpha = .6,
            shape = 21
        ) +
        # geom_point(aes_string(x=x, y=y, color='color'), shape=1) +
        # scale_color_manual(values=color) +
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
        # theme(legend.position = "none") +
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
    sregion <- head(vregion, length(groups))
    sregion$count = singletons[sregion$name, "count"]
    sregion <- sregion %>% mutate(name = paste0(name, " singletons"))
    vregion <- vregion %>% mutate(
        count_perc = round(count / sum(count) * 100, 1),
        count_str = paste0(count, " (", count_perc, "%)")
    )

    # Align the catagory labels
    cat_nudge_y <- 0
    if (length(groups) == 3) { cat_nudge_y <- c(-400, 0, -400) }
    # Shift Count labels
    count_nudge_y <- -10
    if (length(groups) == 3) { count_nudge_y <- c(20, -20, 20, rep(0, nrow(vregion) - 3))  }
    # Shift the singletons stat labels
    label_nudge_y <- 60
    if (length(groups) == 3) { label_nudge_y <- c(60, -60, -60) }

    venn_p <- ggplot() +
        # 1. region count layer
        geom_sf(aes(fill = count), data = venn_region(vdata)) +
        # 2. set edge layer
        # geom_sf(aes(color = factor(id)), data = venn_setedge(data), show.legend = FALSE) +
        # 3. set label layer
        geom_sf_text(aes(label = name), data = venn_setlabel(vdata), nudge_y = cat_nudge_y) +
        # 4. region label layer
        geom_sf_label(
            aes(label = count_str),
            alpha = .8,
            label.padding = unit(.2, "lines"),
            data = vregion,
            nudge_y = count_nudge_y
        ) +
        # 5. singletons label layer
        scale_fill_distiller(palette = "Oranges", direction = 1) +
        new_scale_fill() +
        geom_sf_label(
            aes(label = count, fill = name),
            alpha = .6,
            data = sregion,
            nudge_y = label_nudge_y,
            label.padding = unit(1, "lines"),
            label.r = unit(1.2, "lines"),
            label.size = 0.05,
            show.legend = TRUE
        ) +
        theme_void() +
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        scale_fill_brewer(palette = "Reds", name = "Singletons")

    venn_p
}

plot_upset <- function(counts, singletons) {
    query_singleton <- function(row) { row["Singletons"] == "true" }

    cnts <- column_to_rownames(counts, "CDR3.aa") %>%
        mutate(across(everything(), ~ as.integer(as.logical(.x))))
    sgltns <- unlist(singletons$CDR3.aa)
    cnts$Singletons <- "none"
    cnts[sgltns, "Singletons"] <- "true"
    sets <- setdiff(colnames(cnts), "Singletons")

    upset(cnts, sets = sets, query.legend = "top", sets.x.label = "# clones", queries = list(
        list(
            query = query_singleton,
            color = "orange",
            active = TRUE,
            query.name = "Singletons"
        )
    ))
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
    casedir = file.path(outdir, slugify(casename, tolower = FALSE))
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

    log_info("Handling {i} {case$subject} ...")

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
        log_warn("{casename}, Subject doesn't have enough groups: {subject}")
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

    # # Save samples to group_by so they can be aligned accordingly in the report
    # if (!is.null(section)) {
    #     group_dir <- file.path(casedir, "section")
    #     dir.create(group_dir, showWarnings = FALSE)

    #     sgroups <- subject_row %>%
    #         left_join(cldata) %>%
    #         pull(section) %>%
    #         unique() %>%
    #         paste(collapse = "-")
    #     group_file <- file.path(group_dir, paste0(slugify(sgroups), ".txt"))
    #     cat(subject, file = group_file, sep = "\n", append = TRUE)
    # }

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
            log_warn(
                paste0(
                    "- One of the comparisons doesn't exist in case (", casename,
                    ") for subject (", subject, "): ",
                    case$order[j]
                )
            )
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

        add_report(
            list(
                name = paste0(subject, " (", pair[1], " - ", pair[2], ")"),
                src = scatter_png
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
    png(venn_png, res = 100, height = 600, width = 800)
    print(plot_venndg(counts, groups, singletons))
    dev.off()

    h <- headings(case$section, casename, "Overlapping Clones (Venn Diagram)")
    add_report(
        list(src = venn_png),
        h1 = h$h1,
        h2 = h$h2,
        h3 = h$h3,
        ui = "table_of_images:3"
    )

    upset_dir <- file.path(casedir, "upset")
    upset_png <- file.path(upset_dir, paste0("upset_", slugify(subject), ".png"))
    png(upset_png, res = 100, height = 600, width = 800)
    print(plot_upset(counts, singletons))
    dev.off()

    h <- headings(case$section, casename, "Overlapping Clones (UpSet Plots)")
    add_report(
        list(src = upset_png),
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
