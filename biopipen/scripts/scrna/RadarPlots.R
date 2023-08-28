library(Seurat)
library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggradar)

# input/output
srtfile = {{in.srtobj | r}}
outdir = {{out.outdir | r}}

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

# Used for saving sections
sections = list()

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

run_one_case = function(casename) {
    case = newcases[[casename]]
    print(paste("- Running for case:", casename))

    # Save the section
    if (is.character(case$section) && nchar(case$section) > 0) {
        sections[[case$section]] <<- c(sections[[case$section]], casename)
    }

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

    # Save the counts
    write.table(
        counts,
        file.path(outdir, paste0(casename, ".counts.tsv")),
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
    write.table(
        counts,
        file.path(outdir, paste0(casename, ".percentages.tsv")),
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
    p = ggradar(
        counts %>% as.data.frame() %>% rownames_to_column("group"),
        values.radar = paste0(breaks, "%"),
        grid.min = breaks[1] / 100,
        grid.mid = breaks[2] / 100,
        grid.max = breaks[3] / 100,
        plot.title = casename
    )
    png(
        file.path(outdir, paste0(casename, ".png")),
        width = case$devpars$width,
        height = case$devpars$height,
        res = case$devpars$res
    )
    print(p)
    dev.off()
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

print("- Reading srtobj ...")
srtobj = readRDS(srtfile)
meta = srtobj@meta.data

print("- Mutating meta data if needed ...")
if (is.list(mutaters) && length(mutaters) > 0) {
    mutaters = lapply(mutaters, function(x) parse_expr(x))
    meta = meta %>% mutate(!!!mutaters)
}

# Expand the cases
newcases = list()
for (key in names(cases)) {
    if (is.null(cases[[key]]$each)) {
        newcases[[key]] = cases[[key]]
    } else {
        each_values = meta %>% pull(!!sym(cases[[key]]$each)) %>% unique() %>% na.omit()
        for (evalue in each_values) {
            ekey = if (key == DEFAULT_CASE) evalue else paste0(key, "_", evalue)
            newcases[[ekey]] = cases[[key]]
            newcases[[ekey]]$each_value = evalue
            if (!is.null(cases[[key]]$section)) {
                warn(
                    sprintf("Case %s: `section` is ignored when `each` is specified.", key),
                    immediate. = TRUE
                )
            }
            newcases[[ekey]]$section = key
        }
    }
}

casenames = names(newcases)
sapply(casenames, run_one_case)

print("- Saving sections if any ...")
if (length(sections) > 0) {
    # Write as TOML
    #  section1 = ["case1", "case2"]
    #  section2 = ["case3", "case4"]
    #  ...
    outstr = c()
    for (sec in names(sections)) {
        sec_str = paste0(sec, " = ")
        sec_str = paste0(sec_str, "['", paste(sections[[sec]], collapse = "', '"), "']")
        outstr = c(outstr, sec_str)
    }
    section_file = file.path(outdir, "sections.toml")
    writeLines(outstr, section_file)
}