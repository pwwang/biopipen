source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/mutate_helpers.R")
library(Seurat)
library(rlang)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggVennDiagram)
library(UpSetR)

srtfile <- {{in.srtobj | r}}  # nolint
outdir <- {{out.outdir | r}}  # nolint
mutaters <- {{envs.mutaters | r}}  # nolint
group_by <- {{envs.group_by | r}}  # nolint
group_order <- {{envs.group_order | r}}  # nolint
cells_by <- {{envs.cells_by | r}}  # nolint
cells_order <- {{envs.cells_order | r}}  # nolint
cells_orderby <- {{envs.cells_orderby | r}}  # nolint
cells_n <- {{envs.cells_n | r}}  # nolint
subset <- {{envs.subset | r}}  # nolint
devpars <- {{envs.devpars | r}}  # nolint
each <- {{envs.each | r}}  # nolint
section <- {{envs.section | r}}  # nolint
overlap <- {{envs.overlap | r}}  # nolint
cases <- {{envs.cases | r}}  # nolint

if (is.null(overlap)) { overlap = c() }
overlaps <- list()
print("- Loading seurat object ...")
srtobj <- readRDS(srtfile)

if (!is.null(mutaters) && length(mutaters) > 0) {
    print("- Mutating seurat object ...")
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
                section = section,
                subset = subset
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
                section = section,
                subset = subset
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
            for (ea in eachs) {
                by <- make.names(paste0(".", name, "_", case$each,"_", ea))
                srtobj@meta.data <<- srtobj@meta.data %>%
                    mutate(!!sym(by) := if_else(
                        !!sym(case$each) == ea,
                        !!sym(case$group_by),
                        NA
                    ))
                key <- paste0(case$each, ":", name, " (", ea, ")")
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
    cells_by <- trimws(strsplit(case$cells_by, ",")[[1]])

    sec_case_names <- strsplit(name, ":")[[1]]
    sec_dir <- file.path(outdir, sec_case_names[1])
    casename <- paste(sec_case_names[-1], collapse = ":")
    dir.create(sec_dir, showWarnings = FALSE, recursive = TRUE)
    outfile <- file.path(sec_dir, paste0("case-", casename, ".png"))
    txtfile <- file.path(sec_dir, paste0("case-", casename, ".txt"))

    # subset the seurat object
    meta <- srtobj@meta.data
    if (!is.null(case$subset) && nchar(case$subset) > 0) {
        meta <- dplyr::filter(meta, !!!parse_exprs(case$subset))
    }
    meta <- meta %>%
        drop_na(case$group_by) %>%
        dplyr::filter(!if_all(all_of(cells_by), is.na))

    if (nrow(meta) == 0) {
        stop(paste0("No cells left after filtering NAs for group_by and cells_by for case: ", name))
    }

    if (length(cells_by) > 1) {
        new_cells_by <- paste0(".", paste(cells_by, collapse = "_"))
        meta1 <- meta %>% drop_na(cells_by[1])
        meta1[[new_cells_by]] <- meta1[[cells_by[1]]]
        for (i in 2:length(cells_by)) {
            meta2 <- meta %>% drop_na(cells_by[i])
            meta2[[new_cells_by]] <- meta2[[cells_by[i]]]
            meta1 <- rbind(meta1, meta2)
        }

        cells_by <- new_cells_by
        meta <- meta1
    }

    if (sec_case_names[1] %in% overlap) {
        if (is.null(overlaps[[sec_case_names[1]]])) {
            overlaps[[sec_case_names[1]]] <<- list()
        }
        overlaps[[sec_case_names[1]]][[casename]] <<- meta %>% pull(case$cells_by) %>% unique()
    }

    # add sizes
    meta <- meta %>%
        add_count(!!sym(cells_by), name = "CloneSize") %>%
        add_count(!!sym(cells_by), !!sym(case$group_by), name = "CloneGroupSize") %>%
        add_count(!!sym(cells_by), seurat_clusters, name = "CloneClusterSize") %>%
        add_count(!!sym(cells_by), !!sym(case$group_by), seurat_clusters, name = "CloneGroupClusterSize")

    # filter group_by values not in group_order
    if (!is.null(case$group_order) && length(case$group_order) > 0) {
        meta <- meta %>%
            dplyr::filter(!!sym(case$group_by) %in% case$group_order) %>%
            mutate(!!sym(case$group_by) := factor(!!sym(case$group_by), levels = case$group_order)) %>%
            arrange(!!sym(case$group_by))

        if (nrow(meta) == 0) {
            stop(paste0(
                "No items in `group_order` (",
                paste0(case$group_order, collapse=", "),
                ") in column `", case$group_by , "` for case: ",
                name,
                ". Did you specify the correct `group_by` and `group_order`?"
            ))
        }
    } else {
        meta <- meta %>% arrange(!!sym(case$group_by))
    }

    if (!is.null(case$cells_order) && length(case$cells_order) > 0) {
        # filter cells_by values not in cells_order
        meta <- meta %>%
            dplyr::filter(!!sym(cells_by) %in% case$cells_order) %>%
            mutate(!!sym(cells_by) := factor(!!sym(cells_by), levels = case$cells_order))
    } else if (!is.null(case$cells_orderby)) {
        # otherwise use cells_orderby to order cells_by
        ordered_meta <- meta %>%
            arrange(!!!parse_exprs(case$cells_orderby))
        cells <- ordered_meta %>% pull(cells_by) %>% unique() %>% head(case$cells_n)
        meta <- ordered_meta %>% dplyr::filter(!!sym(cells_by) %in% cells)
        meta[[cells_by]] = factor(meta[[cells_by]], levels = cells)
    }

    write.table(
        meta,
        txtfile,
        sep = "\t",
        row.names = TRUE,
        col.names = TRUE,
        quote = FALSE
    )

    nrows <- length(unique(meta[[cells_by]]))
    ncols <- length(unique(meta[[case$group_by]]))

    devpars <- case$devpars
    if (is.null(devpars)) { devpars = list() }
    if (is.null(devpars$res)) { devpars$res = 100 }
    if (is.null(devpars$width)) { devpars$width = ncols * 100 + 240 }
    if (is.null(devpars$height)) { devpars$height = max(nrows * 100, 600) }

    # plot
    p = meta %>% ggplot(
            aes(
                x = sqrt(CloneGroupSize)/2,
                y = CloneSize,
                width = sqrt(CloneGroupSize),
                fill = seurat_clusters
            )
        ) +
        geom_col(width=.01, position="fill", color = "#888888") +
        geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
        coord_polar("y", start = 0) +
        scale_fill_ucscgb(name = "Cluster", alpha = 1, limits = levels(all_clusters)) +
        theme_void() +
        theme(
            plot.margin = unit(c(1,1,1,1), "cm"),
            legend.text = element_text(size=8),
            legend.title = element_text(size=10)
        )  +
        facet_grid(vars(!!sym(cells_by)), vars(!!sym(case$group_by)), switch="y")

    png(outfile, res = devpars$res, width = devpars$width, height = devpars$height)
    print(p)
    dev.off()
}

do_overlap <- function(section) {
    print(paste("- Running overlaps for section:", section))
    overlap_cases <- overlaps[[section]]
    if (length(overlap_cases) < 2) {
        stop(paste0("Not enough cases for overlap for section: ", section))
    }

    sec_dir <- file.path(outdir, section)
    venn_plot <- file.path(sec_dir, "venn.png")
    venn_p <- ggVennDiagram(overlap_cases, label_percent_digit = 1) +
        scale_fill_distiller(palette = "Reds", direction = 1) +
        scale_x_continuous(expand = expansion(mult = .2))
    png(venn_plot, res = 100, width = 1000, height = 600)
    print(venn_p)
    dev.off()

    upset_plot <- file.path(sec_dir, "upset.png")
    upset_p <- upset(fromList(overlap_cases))
    png(upset_plot, res = 100, width = 800, height = 600)
    print(upset_p)
    dev.off()
}

cases <- expand_cases()
sapply(sort(names(cases)), function(name) do_case(name, cases[[name]]))
sapply(sort(names(overlaps)), do_overlap)
