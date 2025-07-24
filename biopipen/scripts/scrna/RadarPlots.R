{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "repr.R" | source_r }}

library(Seurat)
library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggradar)
library(ggprism)
library(glue)
library(gglogger)

# input/output
srtfile = {{in.srtobj | r}}
outdir = {{out.outdir | r}}
joboutdir = {{job.outdir | r}}

# envs
mutaters = {{envs.mutaters | r}}
by = {{envs.by | r}}
each = {{envs.each | r}}
prefix_each = {{envs.prefix_each | r}}
order = {{envs.order | r}}
colors = {{envs.colors | r}}
ident = {{envs.ident | r}}
cluster_order = {{envs.cluster_order | r}}
breaks = {{envs.breaks | r}}
breakdown = {{envs.breakdown | r}}
test = {{envs.test | r}}
direction = {{envs.direction | r}}
section = {{envs.section | r}}
subset_ = {{envs.subset | r}}
bar_devpars = {{envs.bar_devpars | r}}
devpars = {{envs.devpars | r}}
cases = {{envs.cases | r}}

# DEFAULT_CASE = "DEFAULT"
# sections = c()

log_info("- Reading srtobj ...")
srtobj = biopipen.utils::read_obj(srtfile)
meta = srtobj@meta.data

log_info("- Mutating meta data if needed ...")
if (is.list(mutaters) && length(mutaters) > 0) {
    mutaters <- lapply(mutaters, function(x) parse_expr(x))
    meta <- meta %>% mutate(!!!mutaters)
}

defaults <- list(
    by = by,
    each = each,
    prefix_each = prefix_each,
    order = order,
    colors = colors,
    ident = ident,
    cluster_order = cluster_order,
    breaks = breaks,
    breakdown = breakdown,
    test = test,
    direction = direction,
    section = section,
    subset = subset_,
    bar_devpars = bar_devpars,
    devpars = devpars
)

expand_each <- function(name,  case) {
    outcases <- list()
    if (is.null(case$each) || nchar(case$each) == 0) {
        if (is.null(case$section) || case$section == "DEFAULT") {
            outcases[[name]] <- case
        } else {
            outcases[[paste0(case$section, "::", name)]] <- case
        }
    } else {
        if (is.null(case$subset)) {
            eachs <- meta %>%
                pull(case$each) %>% unique() %>% na.omit() %>% as.vector()
        } else {
            eachs <- meta %>% filter(!!parse_expr(case$subset)) %>%
                pull(case$each) %>% unique() %>% na.omit() %>% as.vector()
        }
        for (each in eachs) {
            if (isTRUE(case$prefix_each)) {
                key <- paste0(name, "::", case$each, " - ", each)
            } else {
                key <- paste0(name, "::", each)
            }
            outcases[[key]] <- case
            outcases[[key]]$section <- name
            outcases[[key]]$each_value <- each
        }
    }
    outcases
}

log_info("- Expanding cases ...")
cases <- expand_cases(cases, defaults, expand_each)

auto_breaks = function(maxval) {
    if (maxval <= 0.1) {  # 10%
        c(0, 5, 10)
    } else if (maxval <= 0.2) {
        c(0, 10, 20)
    } else if (maxval <= 0.3) {
        c(0, 15, 30)
    } else if (maxval <= 0.4) {
        c(0, 20, 40)
    } else if (maxval <= 0.5) {
        c(0, 25, 50)
    } else if (maxval <= 0.6) {
        c(0, 30, 60)
    } else if (maxval <= 0.7) {
        c(0, 35, 70)
    } else if (maxval <= 0.8) {
        c(0, 40, 80)
    } else if (maxval <= 0.9) {
        c(0, 45, 90)
    } else {
        c(0, 50, 100)
    }
}

do_radarplot <- function(info, case, counts) {
    rdr_data = counts %>%
        group_by(!!sym(case$ident), !!sym(case$by)) %>%
        count() %>%
        pivot_wider(
            id_cols = case$by,
            names_from = !!sym(case$ident),
            values_from = n,
            values_fill = 0
        ) %>%
        column_to_rownames(case$by)
        # [by] <cluster1> <cluster2> ...
        #  A       10         20   ...
        #  B       20         30   ...
        #  ...

    # Reorder the clusters if needed
    if (!is.null(case$cluster_order) && length(case$cluster_order) > 0) {
        rdr_data = rdr_data[, case$cluster_order]
    }

    # If clusters are numbers, add a prefix "Cluster"
    if (all(grepl("^\\d+$", colnames(rdr_data)))) {
        colnames(rdr_data) = paste0("Cluster", colnames(rdr_data))
    }

    if (!is.null(case$order) && length(case$order) > 0) {
        rdr_data = rdr_data[case$order, ]
        if (nrow(rdr_data) == 0) {
            stop("No data after reordering. Are items in `order` correct?")
        }
    }

    # Save the counts
    counts_file = file.path(info$casedir, "counts.tsv")
    write.table(
        t(rdr_data),
        counts_file,
        sep = "\t",
        quote = FALSE,
        col.names = TRUE,
        row.names = TRUE
    )

    # Calculate the percentage
    rdr_data = as.matrix(rdr_data)
    if (case$direction == "inter-cluster") {
        rdr_data = t(t(rdr_data) / rowSums(t(rdr_data)))
    } else {
        rdr_data = rdr_data / rowSums(rdr_data)
    }

    # Save the percentages
    perc_file = file.path(info$casedir, "percentages.tsv")
    write.table(
        t(rdr_data),
        perc_file,
        sep = "\t",
        quote = FALSE,
        col.names = TRUE,
        row.names = TRUE
    )

    # Get the breaks
    breaks = if (is.null(case$breaks) || length(case$breaks) == 0) {
        auto_breaks(max(rdr_data))
    } else {
        case$breaks
    }

    # Plot
    if (!is.null(case$colors) && length(case$colors) == 1 && case$colors == "biopipen") {
        colors = pal_biopipen()(nrow(rdr_data))
    } else if (!is.null(case$colors) && length(case$colors) > 0) {
        colors = trimws(unlist(strsplit(case$colors, ",")))
    }

    plotdf <- rdr_data %>%
        as.data.frame() %>%
        rownames_to_column("group") %>%
        mutate(group = factor(group, levels = rownames(rdr_data)))

    p = ggradar(
        plotdf,
        values.radar = paste0(breaks, "%"),
        grid.min = breaks[1] / 100,
        grid.mid = breaks[2] / 100,
        grid.max = breaks[3] / 100,
        group.colours = colors
    )
    prefix <- file.path(info$casedir, "plot")
    save_plot(p, prefix, case$devpars)

    code_file <- paste0(prefix, ".R")
    code = glue(
        "library(ggradar)

        plotdf <- {repr(plotdf)}
        breaks <- {repr(breaks)}
        colors <- {repr(colors)}

        ggradar(
            plotdf,
            values.radar = paste0(breaks, '%'),
            grid.min = breaks[1] / 100,
            grid.mid = breaks[2] / 100,
            grid.max = breaks[3] / 100,
            group.colours = colors
        )"
    )
    writeLines(code, code_file)
}

do_barplot_and_tests <- function(info, case, counts) {
    bardata <- counts %>%
        group_by(!!sym(case$by), !!sym(case$breakdown), !!sym(case$ident)) %>%
        summarise(.n = n(), .groups = "drop")

    write.table(
        bardata,
        file.path(info$casedir, "breakdown-counts.txt"),
        sep = "\t",
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE
    )
    if (case$direction == "inter-cluster") {
        bardata <- bardata %>%
            group_by(!!sym(case$ident)) %>%
            mutate(.frac = .n / sum(.n)) %>%
            ungroup()
    } else {
        bardata <- bardata %>%
            group_by(!!sym(case$by), !!sym(case$breakdown)) %>%
            mutate(.frac = .n / sum(.n)) %>%
            ungroup()
    }

    # Save the percentages
    write.table(
        bardata,
        file.path(info$casedir, "breakdown-percentages.txt"),
        sep = "\t",
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE
    )

    # Reorder the clusters if needed
    if (!is.null(case$cluster_order) && length(case$cluster_order) > 0) {
        bardata <- bardata %>%
            mutate(!!sym(case$ident) := factor(!!sym(case$ident), levels = case$cluster_order))
    }
    # Calculate the mean, mean-sd, mean+sd
    plotdata <- bardata %>%
        group_by(!!sym(case$by), !!sym(case$ident)) %>%
        summarise(.mean = mean(.frac), .sd = sd(.frac), .groups = "drop") %>%
        rowwise() %>%
        mutate(mean_sd1 = max(.mean - .sd, 0), mean_sd2 = .mean + .sd)

    if (!is.null(case$colors) && length(case$colors) == 1 && case$colors == "biopipen") {
        colors <- pal_biopipen(.8)(length(unique(plotdata[[case$by]])))
    } else if (!is.null(case$colors) && length(case$colors) > 0) {
        colors <- trimws(unlist(strsplit(case$colors, ",")))
    }

    # Plot the barplot
    p = ggplot(plotdata, aes(x = !!sym(case$ident), y = .mean, fill = !!sym(case$by))) +
        geom_bar(stat = "identity", position = "dodge", color = "#333333") +
        geom_errorbar(
            aes(ymin = mean_sd1, ymax = mean_sd2),
            width = 0.2,
            alpha = 0.5,
            linewidth = 0.6,
            position = position_dodge(0.9),
            color = "#333333"
        ) +
        theme_prism(axis_text_angle = 45) +
        ylab("Fraction of cells") +
        scale_fill_manual(values = colors)

    prefix = file.path(info$casedir, "barplot")
    save_plot(p, prefix, case$bar_devpars)
    neat_case <- list(by = case$by, ident = case$ident)
    save_plotcode(
        p,
        setup = c(
            'library(rlang)',
            'library(ggplot2)',
            'library(ggprism)',
            '',
            'load("data.RData")',
            'case <- neat_case'
        ),
        prefix,
        "plotdata", "neat_case", "colors")

    # Do the tests in each cluster between groups on .frac
    bys <- bardata %>% pull(!!sym(case$by)) %>% unique()
    if (!is.null(case$test) && test != "none") {
        if (length(bys) < 2) {
            stop("  Cannot do tests with only one group.")
        }

        pairs <- combn(bys, 2, simplify = FALSE)
        test_results <- NULL
        for (pair in pairs) {
            dat <- bardata %>%
                filter(!!sym(case$by) %in% pair) %>%
                select(!!sym(case$by), !!sym(case$ident), .frac) %>%
                group_by(!!sym(case$ident)) %>%
                summarise(
                    comparison = paste0(pair, collapse = " - "),
                    n = paste(as.list(table(!!sym(case$by)))[pair], collapse = "; "),
                    mean = paste(
                        (tibble(.frac, !!sym(case$by)) %>%
                            group_by(!!sym(case$by)) %>%
                            summarise(mean = mean(.frac)) %>%
                            column_to_rownames(case$by) %>%
                            t() %>%
                            as.data.frame() %>%
                            as.list())[pair] %>% unlist() %>% round(3),
                        collapse = "; "
                    ),
                    !!sym(paste0(case$test, "_pval")) := ifelse(
                        case$test == "wilcox",
                        tryCatch(wilcox.test(.frac ~ !!sym(case$by))$p.value, error = function(e) NA),
                        tryCatch(t.test(.frac ~ !!sym(case$by))$p.value, error = function(e) NA)
                    )
                )
            test_results <- rbind(test_results, dat)
        }
        write.table(
            test_results,
            file.path(info$casedir, "tests.txt"),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE
        )
    }
}

add_case_report = function(info, breakdown, test) {
    report = list(
        list(
            name = "Radar Plot",
            contents = list(
                list(
                    kind = "image",
                    src = file.path(info$casedir, "plot.png"),
                    download = list(
                        file.path(info$casedir, "plot.pdf"),
                        list(
                            src = file.path(info$casedir, "plot.R"),
                            tip = "Download the code used to reproduce the plot",
                            icon = "Code"))
                )
            )
        ),
        list(
            name = "Count Table",
            contents = list(
                list(
                    kind = "table",
                    data = list(index_col = 0),
                    src = file.path(info$casedir, "counts.tsv")
                )
            )
        ),
        list(
            name = "Percentage Table",
            contents = list(
                list(
                    kind = "table",
                    data = list(index_col = 0),
                    src = file.path(info$casedir, "percentages.tsv")
                )
            )
        )
    )
    if (!is.null(breakdown)) {
        report = c(
            report,
            list(list(
                name = "Barplot",
                contents = list(
                    list(
                        kind = "image",
                        src = file.path(info$casedir, "barplot.png"),
                        download = list(
                            file.path(info$casedir, "barplot.pdf"),
                            list(
                                src = file.path(info$casedir, "barplot.code.zip"),
                                tip = "Download the code used to reproduce the plot",
                                icon = "Code"
                            )
                        )
                    )
                )
            ))
        )
        if (!is.null(test) && test != "none") {
            report = c(
                report,
                list(list(
                    name = "Tests",
                    contents = list(
                        list(
                            kind = "table",
                            src = file.path(info$casedir, "tests.txt")
                        )
                    )
                ))
            )
        }
    }
    report$h1 = info$h1
    report$h2 = info$h2
    report$ui = "tabs"
    do_call(add_report, report)
}

run_one_case <- function(casename) {
    info <- casename_info(casename, cases, outdir, create = TRUE)
    case <- cases[[casename]]
    log_info("- Running for case: {casename}")

    if (!is.null(case$subset)) {
        m <- meta %>% dplyr::filter(!!rlang::parse_expr(case$subset))
    } else {
        m <- meta
    }
    # Get the counts
    if (!is.null(case$each)) {
        counts <- m %>% dplyr::filter(!!sym(case$each) == case$each_value)
    } else {
        counts <- m
    }
    counts <- counts %>% drop_na(!!sym(case$by)) %>% drop_na(!!sym(case$ident))
    do_radarplot(info, case, counts)

    if (!is.null(case$breakdown)) {
        do_barplot_and_tests(info, case, counts)
    }

    add_case_report(info, case$breakdown, case$test)
}

sapply(sort(names(cases)), run_one_case)

save_report(joboutdir)
