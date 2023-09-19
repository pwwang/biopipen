source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/plot.R")
source("{{biopipen_dir}}/utils/mutate_helpers.R")
library(Seurat)
library(rlang)
library(tidyr)
library(dplyr)
library(ggsci)

srtfile <- {{in.srtobj | r}}  # nolint
outdir <- {{out.outdir | r}}  # nolint
mutaters <- {{envs.mutaters | r}}  # nolint
group_by <- {{envs.group_by | r}}  # nolint
group_order <- {{envs.group_order | r}}  # nolint
cells_by <- {{envs.cells_by | r}}  # nolint
cells_order <- {{envs.cells_order | r}}  # nolint
cells_orderby <- {{envs.cells_orderby | r}}  # nolint
cells_n <- {{envs.cells_n | r}}  # nolint
devpars <- {{envs.devpars | r}}  # nolint
each <- {{envs.each | r}}  # nolint
section <- {{envs.section | r}}  # nolint
cases <- {{envs.cases | r}}  # nolint

srtobj <- readRDS(srtfile)

if (!is.null(mutaters) && length(mutaters) > 0) {
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

all_clusters = srtobj@meta.data %>% pull(seurat_clusters)
if (!is.factor(all_clusters)) {
    all_clusters = factor(all_clusters, levels = sort(unique(all_clusters)))
}

expand_cases <- function() {
    # fill up cases with missing parameters
    if (is.null(cases) || length(cases) == 0) {
        filled_cases <- list(
            DEFAULT = list(
                group_by = group_by,
                group_order = group_order,
                cells_by = cells_by,
                cells_order = cells_order,
                cells_orderby = cells_orderby,
                cells_n = cells_n,
                devpars = devpars,
                each = each,
                section = section
            )
        )
    } else {
        filled_cases <- list()
        for (name in names(cases)) {
            case <- list_setdefault(
                cases[[name]],
                group_by = group_by,
                group_order = group_order,
                cells_by = cells_by,
                cells_order = cells_order,
                cells_orderby = cells_orderby,
                cells_n = cells_n,
                devpars = devpars,
                each = each,
                section = section
            )
            case$devpars <- list_setdefault(case$devpars, devpars)
            filled_cases[[name]] <- case
        }
    }

    outcases <- list()
    # expand each
    for (name in names(filled_cases)) {
        case <- filled_cases[[name]]
        if (is.null(case$each) || nchar(case$each) == 0) {
            outcases[[paste0(case$section, ":", name)]] <- case
        } else {
            eachs <- srtobj@meta.data %>% pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
            for (each in eachs) {
                by <- make.names(paste0(".", name, "_", case$each,"_", each))
                srtobj@meta.data <<- srtobj@meta.data %>%
                    mutate(!!sym(by) := if_else(
                        !!sym(case$each) == each,
                        !!sym(case$group.by),
                        NA
                    ))
                key <- paste0(case$each, ":", name, " (", each, ")")
                outcases[[key]] <- case
                outcases[[key]]$group_by <- by
            }
        }
    }
    outcases
}

do_case <- function(name, case) {
    print(paste("- Running for case:", name))
    if (is.null(case$group_by) || nchar(case$group_by) == 0) {
        stop(paste0("`group_by` must be specified for case", name))
    }
    if (is.null(case$cells_by) || nchar(case$cells_by) == 0) {
        stop(paste0("`cells_by` must be specified for case", name))
    }

    sec_case_names <- strsplit(name, ":")[[1]]
    sec_dir <- file.path(outdir, sec_case_names[1])
    casename <- paste(sec_case_names[-1], collapse = ":")
    dir.create(sec_dir, showWarnings = FALSE, recursive = TRUE)
    outfile <- file.path(sec_dir, paste0(casename, ".png"))
    txtfile <- file.path(sec_dir, paste0(casename, ".txt"))

    # subset the seurat object
    meta <- srtobj@meta.data %>% dplyr::filter(
        !is.na(!!sym(case$group_by)) &
        !is.na(!!sym(case$cells_by))
    )

    # add sizes
    meta <- meta %>%
        add_count(!!sym(case$cells_by), name = "CloneSize") %>%
        add_count(!!sym(case$cells_by), !!sym(case$group_by), name = "CloneGroupSize") %>%
        add_count(!!sym(case$cells_by), seurat_clusters, name = "CloneClusterSize") %>%
        add_count(!!sym(case$cells_by), !!sym(case$group_by), seurat_clusters, name = "CloneGroupClusterSize")

    # filter group_by values not in group_order
    if (!is.null(case$group_order) && length(case$group_order) > 0) {
        meta <- meta %>%
            dplyr::filter(!!sym(case$group_by) %in% case$group_order) %>%
            mutate(!!sym(case$group_by) := factor(!!sym(case$group_by), levels = case$group_order))
    }

    if (!is.null(case$cells_order) && length(case$cells_order) > 0) {
        # filter cells_by values not in cells_order
        meta <- meta %>%
            dplyr::filter(!!sym(case$cells_by) %in% case$cells_order) %>%
            mutate(!!sym(case$cells_by) := factor(!!sym(case$cells_by), levels = case$cells_order))
    } else if (!is.null(case$cells_orderby)) {
        # otherwise use cells_orderby to order cells_by
        ordered_meta <- meta %>% dplyr::arrange(eval(parse(text = case$cells_orderby)))
        cells <- ordered_meta %>% pull(case$cells_by) %>% unique() %>% head(case$cells_n)
        meta <- ordered_meta %>% dplyr::filter(!!sym(case$cells_by) %in% cells)
        meta[[case$cells_by]] = factor(meta[[case$cells_by]], levels = cells)
    }

    write.table(
        meta,
        txtfile,
        sep = "\t",
        row.names = TRUE,
        col.names = TRUE,
        quote = FALSE
    )

    nrows = length(unique(meta[[case$cells_by]]))
    ncols = length(unique(meta[[case$group_by]]))

    devpars <- case$devpars
    if (is.null(devpars)) { devpars = list() }
    if (is.null(devpars$res)) { devpars$res = 100 }
    if (is.null(devpars$width)) { devpars$width = ncols * 100 + 240 }
    if (is.null(devpars$height)) { devpars$height = max(nrows * 100, 600) }

    # plot
    plotGG(
        meta,
        "bar",
        list(
            mapping = aes(
                x = sqrt(CloneGroupSize)/2,
                y = CloneSize,
                width = sqrt(CloneGroupSize),
                fill = seurat_clusters
            ),
            stat = "identity",
            position = "fill"
        ),
        c(
            'geom_col(aes(x=sqrt(CloneGroupSize), y=CloneSize), width=.01, position="fill", color = "#888888")',
            'coord_polar("y", start=0)',
            paste0('facet_grid(vars(', bQuote(case$cells_by), '), vars(', bQuote(case$group_by), '), switch="y")'),
            # 'scale_fill_manual(name = "Cluster", values = pal)',
            # 26-color palette
            'scale_fill_ucscgb(name = "Cluster", limits = levels(all_clusters))',
            'theme_void()',
            'theme(plot.margin = unit(c(1,1,1,1), "cm"))',
            'theme(legend.text=element_text(size=8), legend.title=element_text(size=10))'
        ),
        devpars,
        outfile
    )
}

cases <- expand_cases()
sapply(sort(names(cases)), function(name) do_case(name, cases[[name]]))
