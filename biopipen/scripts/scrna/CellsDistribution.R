source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/mutate_helpers.R")
library(Seurat)
library(rlang)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(UpSetR)
library(slugify)
library(circlize)
library(ComplexHeatmap)

srtfile <- {{in.srtobj | r}}  # nolint
outdir <- {{out.outdir | r}}  # nolint
joboutdir <- {{job.outdir | r}}  # nolint
mutaters <- {{envs.mutaters | r}}  # nolint
cluster_orderby <- {{envs.cluster_orderby | r}}  # nolint
group_by <- {{envs.group_by | r}}  # nolint
group_order <- {{envs.group_order | r}}  # nolint
cells_by <- {{envs.cells_by | r}}  # nolint
cells_order <- {{envs.cells_order | r}}  # nolint
cells_orderby <- {{envs.cells_orderby | r}}  # nolint
cells_n <- {{envs.cells_n | r}}  # nolint
subset <- {{envs.subset | r}}  # nolint
descr <- {{envs.descr | r}}  # nolint
devpars <- {{envs.devpars | r}}  # nolint
each <- {{envs.each | r}}  # nolint
section <- {{envs.section | r}}  # nolint
overlap <- {{envs.overlap | r}}  # nolint
cases <- {{envs.cases | r}}  # nolint

if (is.null(overlap)) { overlap = c() }
overlaps <- list()
log_info("- Loading seurat object ...")
srtobj <- readRDS(srtfile)

if (!is.null(mutaters) && length(mutaters) > 0) {
    log_info("- Mutating seurat object ...")
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

all_clusters = srtobj@meta.data %>% pull(seurat_clusters)
if (!is.factor(all_clusters)) {
    all_clusters = factor(all_clusters, levels = sort(unique(all_clusters)))
}

single_section <- TRUE
expand_cases <- function() {
    # fill up cases with missing parameters
    if (is.null(cases) || length(cases) == 0) {
        filled_cases <- list(
            DEFAULT = list(
                cluster_orderby = cluster_orderby,
                group_by = group_by,
                group_order = group_order,
                cells_by = cells_by,
                cells_order = cells_order,
                cells_orderby = cells_orderby,
                cells_n = cells_n,
                devpars = devpars,
                each = each,
                section = section,
                subset = subset,
                descr = descr
            )
        )
    } else {
        filled_cases <- list()
        for (name in names(cases)) {
            case <- list_setdefault(
                cases[[name]],
                cluster_orderby = cluster_orderby,
                group_by = group_by,
                group_order = group_order,
                cells_by = cells_by,
                cells_order = cells_order,
                cells_orderby = cells_orderby,
                cells_n = cells_n,
                devpars = devpars,
                each = each,
                section = section,
                subset = subset,
                descr = descr
            )
            case$devpars <- list_setdefault(case$devpars, devpars)
            filled_cases[[name]] <- case
        }
    }

    outcases <- list()
    sections <- c()
    # expand each
    for (name in names(filled_cases)) {
        case <- filled_cases[[name]]
        if (is.null(case$each) || nchar(case$each) == 0) {
            sections <- c(sections, case$section)
            outcases[[paste0(case$section, ":", name)]] <- case
        } else {
            sections <- c(sections, case$each)
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
    single_section <<- length(unique(sections)) == 1
    outcases
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
    out$sec_dir <- file.path(outdir, out$section_slug)
    if (create) {
        dir.create(out$sec_dir, showWarnings = FALSE, recursive = TRUE)
    }
    out
}

plot_heatmap <- function(group, meta, case, info, cluster_order_val) {
    log_info(paste("- Running heatmap for case:", info$casename, "group:", group))
    hmfile <- file.path(info$sec_dir, paste0(info$case_slug, ".", slugify(group, tolower = FALSE), ".heatmap.png"))

    # A matrix: 10 × 8 of type int
    #                   g3	g6	g0	g1	g7	g5	g4	g8
    # CSARDATNNEQFF	    8	32	17	26	7	1	NA	NA
    # CASRQNRGSYNEQFF	2	1	20	16	NA	NA	1	NA
    # CSATSYNEQFF	    2	6	3	7	1	8	6	NA
    hmdata <- meta %>% filter(!!sym(case$group_by) == group) %>%
        select(!!sym(case$cells_by), CloneGroupClusterSize, seurat_clusters) %>%
        distinct(!!sym(case$cells_by), seurat_clusters, .keep_all = TRUE) %>%
        pivot_wider(names_from = seurat_clusters, values_from = CloneGroupClusterSize) %>%
        tibble::column_to_rownames(case$cells_by)

    hmdata[, setdiff(levels(meta$seurat_clusters), colnames(hmdata))] <- NA
    # order
    hmdata <- select(hmdata, all_of(levels(meta$seurat_clusters)))

    row_ha <- rowAnnotation(
        Total = anno_barplot(
            hmdata %>% rowSums(na.rm = T),
            gp = gpar(fill = "lightblue", col = NA),
            width = unit(2, "cm")
        )
    )
    ha <- NULL
    extra_height <- 0
    extra_width <- 0  # legend
    if (!is.null(cluster_order_val)) {
        ha <- list()
        ha[[case$cluster_orderby]] <- cluster_order_val
        if (is.numeric(cluster_order_val)) {
            col_fun <- colorRamp2(
                c(min(cluster_order_val), max(cluster_order_val)),
                c("lightyellow", "red"))
            ha$col <- list()
            ha$col[[case$cluster_orderby]] <- col_fun
        }
        ha <- do_call(HeatmapAnnotation, ha)
        extra_height <- 40
        extra_width <- 120
    }
    hm_devpars <- case$hm_devpars
    if (is.null(hm_devpars$res)) { hm_devpars$res = 100 }
    if (is.null(hm_devpars$width)) { hm_devpars$width = ncol(hmdata) * 30 + 400 + extra_width }
    if (is.null(hm_devpars$height)) { hm_devpars$height = nrow(hmdata) * 30 + 60 + extra_height }

    col_fun <- colorRamp2(c(0, max(hmdata, na.rm = T)), c("lightyellow", "purple"))
    png(hmfile, res = hm_devpars$res, width = hm_devpars$width, height = hm_devpars$height)
    p <- Heatmap(
        hmdata,
        name = "Size",
        col = col_fun,
        na_col = "lightyellow",
        row_names_side = "left",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "white", lwd = 1),
        row_names_max_width = max_text_width(rownames(hmdata)),
        right_annotation = row_ha,
        top_annotation = ha
    )
    print(p)
    dev.off()
    hmfile
}

do_case <- function(name, case) {
    log_info(paste("- Running for case:", name))
    if (is.null(case$group_by) || nchar(case$group_by) == 0) {
        stop(paste0("`group_by` must be specified for case", name))
    }
    if (is.null(case$cells_by) || nchar(case$cells_by) == 0) {
        stop(paste0("`cells_by` must be specified for case", name))
    }
    info <- casename_info(name, create = TRUE)
    cells_by <- trimws(strsplit(case$cells_by, ",")[[1]])

    outfile <- file.path(info$sec_dir, paste0(info$case_slug, ".png"))
    txtfile <- file.path(info$sec_dir, paste0(info$case_slug, ".txt"))

    meta <- srtobj@meta.data
    # order the clusters if cluster_orderby is specified
    cluster_order_val <- NULL
    if (!is.null(case$cluster_orderby) && length(case$cluster_orderby) > 0) {
        cluster_order_df <- meta %>%
            group_by(seurat_clusters) %>%
            summarise(
                !!sym(case$cluster_orderby) := !!parse_expr(case$cluster_orderby),
                .groups = "drop") %>%
            arrange(!!sym(case$cluster_orderby))

        cluster_order_val <- pull(cluster_order_df, case$cluster_orderby)

        meta$seurat_clusters <- factor(
            meta$seurat_clusters,
            levels = cluster_order_df %>% pull(seurat_clusters) %>% as.character()
        )
    }

    # subset the seurat object
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
        case$cells_by <- cells_by
        meta <- meta1
    }

    if (info$section %in% overlap) {
        if (is.null(overlaps[[info$section]])) {
            overlaps[[info$section]] <<- list()
        }
        overlaps[[info$section]][[info$case]] <<- meta %>% pull(cells_by) %>% unique()
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
        meta %>% select(
            !!sym(cells_by),
            !!sym(case$group_by),
            seurat_clusters,
            CloneSize,
            CloneGroupSize,
            CloneClusterSize,
            CloneGroupClusterSize,
        ) %>% distinct(
            !!sym(cells_by),
            !!sym(case$group_by),
            seurat_clusters,
            .keep_all = TRUE
        ),
        txtfile,
        sep = "\t",
        row.names = FALSE,
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
        scale_fill_biopipen(name = "Cluster", limits = levels(all_clusters)) +
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

    # heatmaps
    groups = as.character(unique(meta[[case$group_by]]))
    hmfigs = sapply(groups, plot_heatmap, meta, case, info, cluster_order_val)

    add_report(
        list(
            kind = "descr",
            content = ifelse(
                is.null(case$descr) || nchar(case$descr) == 0,
                paste0(
                    "Distribution for cells in ",
                    "<code>", html_escape(cells_by), "</code>",
                    " for ",
                    "<code>", html_escape(case$group_by), "</code>"
                ),
                case$descr
            )
        ),
        h1 = ifelse(
            info$section == "DEFAULT",
            info$case,
            ifelse(single_section, paste0(info$section, " - ", info$case), info$section)
        ),
        h2 = ifelse(single_section, "#", info$case)
    )

    add_report(
        list(
            name = "Pie Charts",
            contents = list(list(
                kind = "image",
                src = outfile
            ))
        ),
        list(
            name = "Heatmaps",
            ui = "table_of_images",
            contents = lapply(seq_along(groups), function(i) {
                list(descr = groups[i], src = hmfigs[i])
            })
        ),
        list(
            name = "Distribution Table",
            contents = list(list(
                kind = "table",
                data = list(nrows = 100),
                src = txtfile
            ))
        ),
        h1 = ifelse(
            info$section == "DEFAULT",
            info$case,
            ifelse(single_section, paste0(info$section, " - ", info$case), info$section)
        ),
        h2 = ifelse(single_section, "#", info$case),
        ui = "tabs"
    )
}

do_overlap <- function(section) {
    log_info(paste("- Running overlaps for section:", section))
    overlap_cases <- overlaps[[section]]
    if (length(overlap_cases) < 2) {
        stop(paste0("Not enough cases for overlap for section: ", section))
    }

    sec_dir <- file.path(outdir, slugify(section, tolower = FALSE))
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

    add_report(
        list(
            name = "Venn Plot",
            contents = list(list(
                kind = "image",
                src = venn_plot
            ))
        ),
        list(
            name = "UpSet Plot",
            contents = list(list(
                kind = "image",
                src = upset_plot
            ))
        ),
        h1 = "Overlapping Groups",
        h2 = section,
        ui = "tabs"
    )
}

cases <- expand_cases()
sapply(sort(names(cases)), function(name) do_case(name, cases[[name]]))
sapply(sort(names(overlaps)), do_overlap)

save_report(joboutdir)
