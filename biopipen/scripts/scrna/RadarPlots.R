source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggradar)
library(slugify)

# input/output
srtfile = {{in.srtobj | r}}
outdir = {{out.outdir | r}}
joboutdir = {{job.outdir | r}}

# envs
mutaters = {{envs.mutaters | r}}
by = {{envs.by | r}}
each = {{envs.each | r}}
order = {{envs.order | r}}
cluster_col = {{envs.cluster_col | r}}
cluster_order = {{envs.cluster_order | r}}
breaks = {{envs.breaks | r}}
direction = {{envs.direction | r}}
section = {{envs.section | r}}
devpars = {{envs.devpars | r}}
cases = {{envs.cases | r}}

DEFAULT_CASE = "DEFAULT"
sections = c()

log_info("Reading srtobj ...")
srtobj = readRDS(srtfile)
meta = srtobj@meta.data

log_info("Mutating meta data if needed ...")
if (is.list(mutaters) && length(mutaters) > 0) {
    mutaters = lapply(mutaters, function(x) parse_expr(x))
    meta = meta %>% mutate(!!!mutaters)
}

# fill up the cases
if (length(cases) == 0) {
    cases[[DEFAULT_CASE]] = list(
        by = by,
        each = each,
        order = order,
        cluster_col = cluster_col,
        cluster_order = cluster_order,
        breaks = breaks,
        direction = direction,
        section = section,
        devpars = devpars
    )
} else {
    # Use the values given directly under `envs` as default
    for (key in names(cases)) {
        if (is.null(cases[[key]]$by)) {
            cases[[key]]$by = by
        }
        if (is.null(cases[[key]]$each)) {
            cases[[key]]$each = each
        }
        if (is.null(cases[[key]]$order)) {
            cases[[key]]$order = order
        }
        if (is.null(cases[[key]]$cluster_col)) {
            cases[[key]]$cluster_col = cluster_col
        }
        if (is.null(cases[[key]]$cluster_order)) {
            cases[[key]]$cluster_order = cluster_order
        }
        if (is.null(cases[[key]]$breaks)) {
            cases[[key]]$breaks = breaks
        }
        if (is.null(cases[[key]]$direction)) {
            cases[[key]]$direction = direction
        }
        if (is.null(cases[[key]]$section)) {
            cases[[key]]$section = section
        }
        if (is.null(cases[[key]]$devpars)) {
            cases[[key]]$devpars = devpars
        }
        if (is.null(cases[[key]]$devpars$width)) {
            cases[[key]]$devpars$width = devpars$width
        }
        if (is.null(cases[[key]]$devpars$height)) {
            cases[[key]]$devpars$height = devpars$height
        }
        if (is.null(cases[[key]]$devpars$res)) {
            cases[[key]]$devpars$res = devpars$res
        }
    }
}

# Expand the cases
newcases = list()
for (key in names(cases)) {
    if (is.null(cases[[key]]$each)) {
        sections <- c(sections, cases[[key]]$section)
        newcases[[paste0(cases[[key]]$section, ":", key)]] = cases[[key]]
    } else {
        each_values = meta %>% pull(!!sym(cases[[key]]$each)) %>% unique() %>% na.omit()
        sections <- c(sections, key)
        for (evalue in each_values) {
            ekey = paste0(key, ":", evalue)
            newcases[[ekey]] = cases[[key]]
            newcases[[ekey]]$each_value = evalue
            if (!is.null(cases[[key]]$section)) {
                log_warn(
                    sprintf("Case %s: `section` is ignored when `each` is specified.", key)
                )
            }
        }
    }
}

single_section <- length(sections) == 1 && sections[[1]] == "DEFAULT"

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

casename_info <- function(casename, create = FALSE) {
    sec_case_names <- strsplit(casename, ":")[[1]]
    cname <- paste(sec_case_names[-1], collapse = ":")

    out <- list(
        casename = casename,
        section = sec_case_names[1],
        case = cname,
        section_slug = slugify(sec_case_names[1], tolower = FALSE),
        case_slug = slugify(cname, tolower = FALSE)
    )
    out$casedir <- file.path(outdir, out$section_slug, out$case_slug)
    if (create) {
        dir.create(out$casedir, showWarnings = FALSE, recursive = TRUE)
    }
    out
}

run_one_case = function(casename) {
    info = casename_info(casename, create = TRUE)
    case = newcases[[casename]]
    log_info("Running for case: {casename}")

    # Get the counts
    counts = if (!is.null(case$each)) meta %>% filter(!!sym(case$each) == case$each_value) else meta
    counts = counts %>%
        filter(!is.na(!!sym(case$by))) %>%
        group_by(!!sym(case$cluster_col), !!sym(case$by)) %>%
        count() %>%
        pivot_wider(
            id_cols = case$by,
            names_from = !!sym(case$cluster_col),
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
        counts = counts[, case$cluster_order]
    }

    # If clusters are numbers, add a prefix "Cluster"
    if (all(grepl("^\\d+$", colnames(counts)))) {
        colnames(counts) = paste0("Cluster", colnames(counts))
    }

    if (!is.null(case$order) && length(case$order) > 0) {
        counts = counts[case$order, ]
        if (nrow(counts) == 0) {
            stop("No data after reordering. Are items in `order` correct?")
        }
    }

    # Save the counts
    counts_file = file.path(info$casedir, "counts.tsv")
    write.table(
        counts,
        counts_file,
        sep = "\t",
        quote = FALSE,
        col.names = TRUE,
        row.names = TRUE
    )

    # Calculate the percentage
    counts = as.matrix(counts)
    if (case$direction == "inter-cluster") {
        counts = t(t(counts) / rowSums(t(counts)))
    } else {
        counts = counts / rowSums(counts)
    }

    # Save the percentages
    perc_file = file.path(info$casedir, "percentages.tsv")
    write.table(
        counts,
        perc_file,
        sep = "\t",
        quote = FALSE,
        col.names = TRUE,
        row.names = TRUE
    )

    # Get the breaks
    breaks = if (is.null(case$breaks) || length(case$breaks) == 0) {
        auto_breaks(max(counts))
    } else {
        case$breaks
    }

    # Plot
    plotfile = file.path(info$casedir, "plot.png")
    p = ggradar(
        counts %>% as.data.frame() %>% rownames_to_column("group"),
        values.radar = paste0(breaks, "%"),
        grid.min = breaks[1] / 100,
        grid.mid = breaks[2] / 100,
        grid.max = breaks[3] / 100,
        group.colours = pal_biopipen()(nrow(counts))
    )
    png(
        plotfile,
        width = case$devpars$width,
        height = case$devpars$height,
        res = case$devpars$res
    )
    print(p)
    dev.off()

    add_case_report(info)
}

add_case_report = function(info) {
    add_report(
        list(
            name = "Radar Plot",
            contents = list(
                list(
                    kind = "image",
                    src = file.path(info$casedir, "plot.png")
                )
            )
        ),
        list(
            name = "Count Table",
            contents = list(
                list(
                    kind = "table",
                    src = file.path(info$casedir, "counts.tsv")
                )
            )
        ),
        list(
            name = "Percentage Table",
            contents = list(
                list(
                    kind = "table",
                    src = file.path(info$casedir, "percentages.tsv")
                )
            )
        ),
        h1 = ifelse(
            info$section == "DEFAULT",
            info$case,
            ifelse(single_section, paste0(info$section, " - ", info$case), info$section)
        ),
        h2 = ifelse(
            info$section == "DEFAULT",
            "#",
            ifelse(single_section, "#", info$case)
        ),
        ui = "tabs"
    )
}


casenames = names(newcases)
sapply(casenames, run_one_case)

save_report(joboutdir)
