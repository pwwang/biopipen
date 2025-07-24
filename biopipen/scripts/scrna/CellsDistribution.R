{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "mutate_helpers.R" | source_r }}

library(Seurat)
library(rlang)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggVennDiagram)
library(UpSetR)
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
hm_devpars <- {{envs.hm_devpars | r}}  # nolint
each <- {{envs.each | r}}  # nolint
section <- {{envs.section | r}}  # nolint
prefix_each <- {{envs.prefix_each | r}}  # nolint
overlap <- {{envs.overlap | r}}  # nolint
cases <- {{envs.cases | r}}  # nolint

overlap <- overlap %||% c()
overlaps <- list()
log_info("- Loading seurat object ...")
srtobj <- biopipen.utils::read_obj(srtfile)

if (!is.null(mutaters) && length(mutaters) > 0) {
    log_info("- Mutating seurat object ...")
    srtobj@meta.data <- srtobj@meta.data %>% mutate(!!!lapply(mutaters, parse_expr))
}

defaults <- list(
    cluster_orderby = cluster_orderby,
    group_by = group_by,
    group_order = group_order,
    cells_by = cells_by,
    cells_order = cells_order,
    cells_orderby = cells_orderby,
    cells_n = cells_n,
    devpars = devpars,
    hm_devpars = hm_devpars,
    each = each,
    section = section,
    prefix_each = prefix_each,
    subset = subset,
    descr = descr
)

expand_each <- function(name, case) {
    outcases <- list()
    if (is.null(case$each) || nchar(case$each) == 0) {
        if (is.null(case$section) || case$section == "DEFAULT") {
            outcases[[name]] <- case
        } else {
            outcases[[paste0(case$section, "::", name)]] <- case
        }
    } else {
        if (!is.null(case$section) && case$section != "DEFAULT") {
            log_warn("Ignoring `section` in case `{name}` when `each` is set.")
            case$section <- NULL
        }
        if (!is.null(case$subset)) {
            eachs <- srtobj@meta.data %>%
                dplyr::filter(!!parse_expr(case$subset)) %>%
                pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
        } else {
            eachs <- srtobj@meta.data %>%
                pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
        }

        for (each in eachs) {
            by <- make.names(paste0(".", name, "_", case$each,"_", each))
            srtobj@meta.data <<- srtobj@meta.data %>%
                mutate(!!sym(by) := if_else(
                    !!sym(case$each) == each,
                    !!sym(case$group_by),
                    NA
                ))
            if (isTRUE(case$prefix_each)) {
                key <- paste0(name, "::", case$each, " - ", each)
            } else {
                key <- paste0(name, "::", each)
            }
            outcases[[key]] <- case
            outcases[[key]]$section <- name
            outcases[[key]]$group_by <- by
        }
    }
    outcases
}

log_info("- Expanding cases...")
cases <- expand_cases(cases, defaults, expand_each)

plot_heatmap <- function(m, cells_by, group_by, cluster_order_val, cluster_orderby) {
    # A matrix: 10 × 8 of type int
    #                   g3	g6	g0	g1	g7	g5	g4	g8
    # CSARDATNNEQFF	    8	32	17	26	7	1	NA	NA
    # CASRQNRGSYNEQFF	2	1	20	16	NA	NA	1	NA
    # CSATSYNEQFF	    2	6	3	7	1	8	6	NA
    hmdata <- m %>%
        mutate(
            !!sym(cells_by) := paste0("[", !!sym(group_by), "] ", !!sym(cells_by))
        ) %>%
        select(!!sym(cells_by), CloneGroupClusterSize, seurat_clusters) %>%
        distinct(!!sym(cells_by), seurat_clusters, .keep_all = TRUE) %>%
        pivot_wider(names_from = seurat_clusters, values_from = CloneGroupClusterSize) %>%
        tibble::column_to_rownames(cells_by)

    hmdata[, setdiff(levels(m$seurat_clusters), colnames(hmdata))] <- NA
    # order
    hmdata <- select(hmdata, all_of(levels(m$seurat_clusters)))

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
        ha[[cluster_orderby]] <- cluster_order_val
        if (is.numeric(cluster_order_val)) {
            col_fun <- colorRamp2(
                c(min(cluster_order_val), max(cluster_order_val)),
                c("lightyellow", "red"))
            ha$col <- list()
            ha$col[[cluster_orderby]] <- col_fun
        }
        ha <- do_call(HeatmapAnnotation, ha)
        extra_height <- 40
        extra_width <- 120
    }

    col_fun <- colorRamp2(c(0, max(hmdata, na.rm = T)), c("lightyellow", "purple"))
    Heatmap(
        as.matrix(hmdata),
        name = cells_by,
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
}

sections <- c()
do_case <- function(name, case) {
    log_info("- Handling case: {name}")
    if (is.null(case$group_by) || nchar(case$group_by) == 0) {
        stop(paste0("  `group_by` must be specified for case", name))
    }
    if (is.null(case$cells_by) || nchar(case$cells_by) == 0) {
        stop(paste0("  `cells_by` must be specified for case", name))
    }

    info <- casename_info(name, cases, outdir, create = TRUE)
    cells_by <- trimws(strsplit(case$cells_by, ",")[[1]])

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
    if (!is.factor(meta$seurat_clusters)) {
        meta$seurat_clusters <- factor(meta$seurat_clusters, levels = sort(unique(meta$seurat_clusters)))
    }
    all_clusters <- meta$seurat_clusters

    # subset the seurat object
    if (!is.null(case$subset) && nchar(case$subset) > 0) {
        meta <- dplyr::filter(meta, !!!parse_exprs(case$subset))
    }
    meta <- meta %>% drop_na(case$group_by)
        # dplyr::filter(!if_all(all_of(cells_by), is.na))

    if (nrow(meta) == 0) {
        stop(paste0("No cells left after filtering NAs for `group_by`"))
    }

    if (info$section %in% overlap && length(cells_by) > 1) {
        stop(paste0("Overlapping groups can only be done for a single `cells_by`"))
    }

    # filter group_by values not in group_order
    if (!is.null(case$group_order) && length(case$group_order) > 0) {
        meta <- meta %>%
            dplyr::filter(!!sym(case$group_by) %in% case$group_order) %>%
            mutate(!!sym(case$group_by) := factor(!!sym(case$group_by), levels = case$group_order))

        if (nrow(meta) == 0) {
            stop(paste0(
                "No items in `group_order` (",
                paste0(case$group_order, collapse=", "),
                ") in column `", case$group_by ,
                ". Did you specify the correct `group_by` and `group_order`?"
            ))
        }
    } else {
        meta <- meta %>% mutate(!!sym(case$group_by) := factor(!!sym(case$group_by)))
    }
    ngroups <- length(unique(meta[[case$group_by]]))
    sections <<- union(sections, info$section)

    piecharts <- list()
    hmdata <- NULL
    hmrowlbls <- c()
    hmsplits <- c()
    hmfile <- file.path(info$casedir, paste0(info$case_slug, ".heatmap.png"))
    cells_rows <- 0
    table_files <- c()
    for (n in seq_along(cells_by)) {
        cby <- cells_by[n]
        log_info("- Processing cells_by: {cby}")
        m <- meta %>% drop_na(!!sym(cby))

        # check if there are enough cells
        if (nrow(m) == 0) {
            stop(paste0("  No cells left after filtering NAs for `", cby, "`"))
        }

        if (info$section %in% overlap) {
            overlaps[[info$section]] <<- overlaps[[info$section]] %||% list()
            overlaps[[info$section]][[info$case]] <<- m %>% pull(cby) %>% unique()
        }

        # add sizes
        m <- m %>%
            add_count(!!sym(cby), name = "CloneSize") %>%
            add_count(!!sym(cby), !!sym(case$group_by), name = "CloneGroupSize") %>%
            add_count(!!sym(cby), seurat_clusters, name = "CloneClusterSize") %>%
            add_count(!!sym(cby), !!sym(case$group_by), seurat_clusters, name = "CloneGroupClusterSize") %>%
            select(
                !!sym(cby),
                !!sym(case$group_by),
                seurat_clusters,
                CloneSize,
                CloneGroupSize,
                CloneClusterSize,
                CloneGroupClusterSize,
            ) %>% distinct(
                !!sym(cby),
                !!sym(case$group_by),
                seurat_clusters,
                .keep_all = TRUE
            )


        # apply cells order
        if (!is.null(case$cells_order) && length(case$cells_order) > 0) {
            m <- m %>%
                dplyr::filter(!!sym(cby) %in% case$cells_order) %>%
                mutate(!!sym(cby) := factor(!!sym(cby), levels = case$cells_order))
        } else if (!is.null(case$cells_orderby)) {
            ordered_m <- m %>% arrange(!!!parse_exprs(case$cells_orderby))
            cells <- ordered_m %>% pull(cby) %>% unique() %>% head(case$cells_n)
            m <- ordered_m %>% dplyr::filter(!!sym(cby) %in% cells)
            m[[cby]] = factor(m[[cby]], levels = cells)
        }

        # save the filtered data
        table_file <- file.path(info$casedir, paste0(info$case_slug, ".", slugify(cby), ".txt"))
        table_files <- c(table_files, table_file)
        write.table(
            m, table_file,
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
        )

        log_debug("  Plotting pie charts ...")
        cells_rows <- cells_rows + length(unique(m[[cby]]))
        if (n == 1) {
            plot.margin <- unit(c(1,1,0,1), "cm")
            strip.text.x <- element_text(margin = margin(b = 0.5, unit = "cm"))
        } else if (n == length(cells_by)) {
            plot.margin <- unit(c(0,1,1,1), "cm")
            strip.text.x <- element_blank()
        } else {
            plot.margin <- unit(c(0,1,0,1), "cm")
            strip.text.x <- element_blank()
        }
        p = m %>% ggplot(
                aes(
                    x = sqrt(CloneGroupSize)/2,
                    y = CloneGroupClusterSize,
                    width = sqrt(CloneGroupSize),
                    fill = seurat_clusters
                )
            ) +
            geom_col(width=.01, position="fill", color = "#888888") +
            geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
            coord_polar("y", start = 0) +
            scale_fill_manual(name = "Cluster", values = pal_biopipen()(length(levels(all_clusters)))) +
            theme_void() +
            theme(
                plot.margin = plot.margin,
                legend.margin = margin(l = .8, unit = "cm"),
                legend.text = element_text(size=6),
                legend.title = element_text(size=8),
                legend.key.size = unit(0.5, "cm"),
                strip.text.x = strip.text.x,
                strip.text.y = element_text(
                    angle = 0, hjust = 1,
                    margin = margin(r = 0.5, unit = "cm"))
            )  +
            facet_grid(vars(!!sym(cby)), vars(!!sym(case$group_by)), switch="y")

        piecharts[[length(piecharts) + 1]] <- p

        # heatmaps
        log_debug("  Preparing pie charts ...")
        hmd <- m %>%
            arrange(!!sym(case$group_by), !!sym(cby)) %>%
            mutate(!!sym(cby) := paste0("[", !!sym(case$group_by), "] ", !!sym(cby))) %>%
            select(!!sym(cby), CloneGroupClusterSize, seurat_clusters) %>%
            distinct(!!sym(cby), seurat_clusters, .keep_all = TRUE) %>%
            pivot_wider(names_from = seurat_clusters, values_from = CloneGroupClusterSize) %>%
            tibble::column_to_rownames(cby)
        hmd[, setdiff(levels(m$seurat_clusters), colnames(hmd))] <- NA
        hmd <- select(hmd, all_of(levels(m$seurat_clusters)))
        hmsplits <- c(hmsplits, rep(cby, nrow(hmd)))
        hmrowlbls <- c(hmrowlbls, rownames(hmd))
        rownames(hmd) <- NULL

        hmdata <- bind_rows(hmdata, hmd)
    }

    log_info("  Merging and saving pie charts ...")
    devpars = case$devpars
    # assemble and save pie chart plots
    res <- devpars$res %||% 100
    #                         legend, cells_by names
    width <- devpars$width %||% (400 + 120 + 100 * ngroups)
    #                         group_by names
    height <- devpars$height %||% (120 + 100 * cells_rows)

    p <- wrap_plots(piecharts, ncol = 1, guides = "collect")

    piefile <- file.path(info$casedir, paste0(info$case_slug, ".png"))
    png(piefile, res = res, width = width, height = height)
    print(p)
    dev.off()

    piefile_pdf <- file.path(info$casedir, paste0(info$case_slug, ".pdf"))
    pdf(piefile_pdf, width = width / res, height = height / res)
    print(p)
    dev.off()

    log_info("  Plotting and saving heatmap ...")
    row_ha <- rowAnnotation(
        Total = anno_barplot(
            hmdata %>% rowSums(na.rm = T),
            gp = gpar(fill = "lightblue", col = NA),
            width = unit(1.5, "cm")
        )
    )
    ha <- NULL
    extra_height <- 0
    extra_width <- 0  # legend
    if (!is.null(cluster_order_val)) {
        ha <- list()
        ha[[cluster_orderby]] <- cluster_order_val
        if (is.numeric(cluster_order_val)) {
            col_fun <- colorRamp2(
                c(min(cluster_order_val), max(cluster_order_val)),
                c("lightyellow", "red"))
            ha$col <- list()
            ha$col[[cluster_orderby]] <- col_fun
        }
        ha <- do_call(HeatmapAnnotation, ha)
        extra_height <- 40
        extra_width <- 120
    }
    if (length(cells_by) == 1) {
        hmsplits <- NULL
        extra_width <- extra_width - 15
    } else {
        # keep the row order
        hmsplits <- factor(hmsplits, levels = unique(hmsplits))
    }

    col_fun <- colorRamp2(c(0, max(hmdata, na.rm = T)), c("lightyellow", "purple"))
    hm_devpars <- case$hm_devpars
    hm_res <- hm_devpars$res %||% 100
    hm_width <- hm_devpars$width %||% (600 + 15 * length(unique(meta$seurat_clusters)) + extra_width)
    hm_height <- hm_devpars$height %||% (450 + 15 * cells_rows + extra_height)
    hm <- Heatmap(
        as.matrix(hmdata),
        name = "Size",
        col = col_fun,
        na_col = "lightyellow",
        row_names_side = "left",
        row_names_max_width = max_text_width(
            hmrowlbls,
            gp = gpar(fontsize = 12)
        ),
        row_labels = hmrowlbls,
        row_split = hmsplits,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "white", lwd = 1),
        right_annotation = row_ha,
        top_annotation = ha
    )
    png(hmfile, res = hm_res, width = hm_width, height = hm_height)
    print(hm)
    dev.off()

    hmfile_pdf <- gsub(".png$", ".pdf", hmfile)
    pdf(hmfile_pdf, width = hm_width / hm_res, height = hm_height / hm_res)
    print(hm)
    dev.off()

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
        h1 = info$h1,
        h2 = info$h2
    )

    add_report(
        list(
            name = "Pie Charts",
            contents = list(list(kind = "image", src = piefile, download = piefile_pdf))
        ),
        list(
            name = "Heatmap",
            contents = list(list(src = hmfile, kind = "image", download = hmfile_pdf))
        ),
        list(
            name = "Distribution Table",
            contents = do.call(c, lapply(
                seq_along(cells_by),
                function(i) list(
                    list(kind = "descr", content = paste0("Cells by: ", cells_by[i])),
                    list(kind = "table", data = list(nrows = 100), src = table_files[i])
                )
            ))
        ),
        h1 = info$h1,
        h2 = info$h2,
        ui = "tabs"
    )
}

do_overlap <- function(section) {
    log_info(paste("- Running overlaps for section:", section))
    overlap_cases <- overlaps[[section]]
    if (length(overlap_cases) < 2) {
        stop(paste0("  Not enough cases for overlap for section: ", section))
    }

    sec_dir <- file.path(outdir, paste0("overlap - ", slugify(section)))
    dir.create(sec_dir, showWarnings = FALSE)
    venn_plot <- file.path(sec_dir, "venn.png")
    venn_p <- ggVennDiagram(overlap_cases, label_percent_digit = 1) +
        scale_fill_distiller(palette = "Reds", direction = 1) +
        scale_x_continuous(expand = expansion(mult = .2))
    png(venn_plot, res = 100, width = 1000, height = 600)
    print(venn_p)
    dev.off()

    venn_plot_pdf <- gsub(".png$", ".pdf", venn_plot)
    pdf(venn_plot_pdf, width = 10, height = 6)
    print(venn_p)
    dev.off()

    upset_plot <- file.path(sec_dir, "upset.png")
    upset_p <- upset(fromList(overlap_cases))
    png(upset_plot, res = 100, width = 800, height = 600)
    print(upset_p)
    dev.off()

    upset_plot_pdf <- gsub(".png$", ".pdf", upset_plot)
    pdf(upset_plot_pdf, width = 8, height = 6)
    print(upset_p)
    dev.off()

    add_report(
        list(
            name = "Venn Plot",
            contents = list(list(
                kind = "image",
                src = venn_plot,
                download = venn_plot_pdf
            ))
        ),
        list(
            name = "UpSet Plot",
            contents = list(list(
                kind = "image",
                src = upset_plot,
                download = upset_plot_pdf
            ))
        ),
        h1 = "Overlapping Groups",
        h2 = section,
        ui = "tabs"
    )
}

sapply(sort(names(cases)), function(name) do_case(name, cases[[name]]))

unhit_overlaps <- setdiff(overlap, names(overlaps))
if (length(unhit_overlaps) > 0) {
    log_warn("- No sections found for overlapping analysis: {paste(unhit_overlaps, collapse=', ')}")
    log_warn("  Available sections: {paste(sections, collapse=', ')}")
}
sapply(sort(names(overlaps)), do_overlap)

save_report(joboutdir)
