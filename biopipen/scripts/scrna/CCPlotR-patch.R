# patched version of cc_circos
# See https://github.com/Sarah145/CCPlotR/issues/4

cc_circos <- function(cc_df, option = "A", n_top_ints = 15, exp_df = NULL, cell_cols = NULL, palette = "BuPu", cex = 1, show_legend = TRUE, scale = FALSE, ...) {
    stopifnot("'cc_df' must be a dataframe" = is(cc_df, "data.frame"))
    stopifnot("cc_df should contain columns named source, target, ligand, receptor and score. See `toy_data` for an example." = all(c('source', 'target', 'ligand', 'receptor', 'score') %in% colnames(cc_df)))
    stopifnot("option must be either 'A', 'B', 'C'" = option %in% c('A', 'B', 'C'))
    library(stringr)
    library(ComplexHeatmap)
    library(circlize)
    circos.clear()

    target <- score <- ligand <- receptor <- source_lig <- target_rec <- cell_type <- gene <- cell_gene <- NULL
    if (option == "A") {
        input_df <- cc_df %>%
            mutate(source = factor(source), target = factor(target)) %>%
            group_by(source, target) %>%
            tally()
        if (is.null(cell_cols)) {
            cell_cols <- setNames(paletteMartin(n = length(unique(c(input_df$source, input_df$target)))), unique(c(input_df$source, input_df$target)))
        }
        circlize_plot <- function() {
            par(cex = cex)
            chordDiagram(input_df,
                scale = FALSE, grid.col = cell_cols,
                annotationTrack = c("grid", "name"), directional = 1, direction.type = c("arrows", "diffHeight"), link.arr.type = "big.arrow", link.arr.length = 0.1, diffHeight = -mm_h(0.5), preAllocateTracks = list(
                    track.height = mm_h(10),
                    track.margin = c(mm_h(2), -mm_h(4))
                ), ...
            )
        }
    } else if (option == "B") {
        input_df <- cc_df %>%
            slice_max(order_by = score, n = n_top_ints) %>%
            mutate(
                source_lig = paste0(source, "|", ligand),
                target_rec = paste0(target, "|", receptor)
            )
        arr_wd <- (((input_df$score - min(input_df$score)) / (max(input_df$score) - min(input_df$score))) * (4)) + 1

        if (is.null(cell_cols)) {
            cell_cols <- setNames(paletteMartin(n = length(unique(c(input_df$source, input_df$target)))), unique(c(input_df$source, input_df$target)))
        }

        link_cols <- c()
        for (i in input_df$source_lig) {
            link_cols <- c(link_cols, cell_cols[str_extract(i, "[^|]+")])
        }

        segments <- unique(c(paste0(input_df$source, "|", input_df$ligand), paste0(input_df$target, "|", input_df$receptor)))
        grp <- str_extract(segments, "[^|]+")
        names(grp) <- segments
        lgd <- Legend(
            labels = unique(c(input_df$source, input_df$target)),
            title = "Cell type",
            type = "points",
            title_gp = gpar(fontsize = 14 * cex),
            labels_gp = gpar(fontsize = 12 * cex),
            legend_gp = gpar(col = "transparent"),
            background = cell_cols[unique(c(input_df$source, input_df$target))]
        )
        circlize_plot <- function() {
            par(cex = cex)
            chordDiagram(
                input_df %>%
                    select(source_lig, target_rec, score),
                directional = 1, group = grp, link.sort = FALSE, scale = scale, diffHeight = 0.005,
                direction.type = c("arrows"), link.arr.type = "triangle", annotationTrack = c(),
                preAllocateTracks = list(list(track.height = 0.175), list(track.height = 0.05)),
                big.gap = 3, transparency = 1, link.arr.lwd = arr_wd, link.arr.col = link_cols,
                link.arr.length = 0.4, link.arr.width = 0.35, ...
            )
            circos.track(track.index = 1, panel.fun = function(x, y) {
                circos.text(CELL_META$xcenter, CELL_META$ylim[1], str_extract(CELL_META$sector.index, "[^|]+$"),
                    facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1.3
                )
            }, bg.border = NA)
            for (l in unique(str_extract(segments, "[^|]+"))) {
                highlight.sector(segments[str_detect(segments, paste0("^", str_escape(l)))], track.index = 2, col = cell_cols[l])
            }
            if (show_legend == TRUE) {
                draw(lgd, just = c("left", "bottom"), x = unit(5, "mm"), y = unit(5, "mm"))
            }
            circos.clear()
        }
    } else if (option == "C") {
        stopifnot("'exp_df' must be a dataframe" = is(exp_df, "data.frame"))
        stopifnot("exp_df should contain columns named cell_type, gene and mean_exp. See `toy_exp` for an example." = all(c('cell_type', 'gene', 'mean_exp') %in% colnames(exp_df)))

        input_df <- cc_df %>%
            slice_max(order_by = score, n = n_top_ints) %>%
            mutate(
                source_lig = paste0(source, "|", ligand),
                target_rec = paste0(target, "|", receptor)
            )

        arr_wd <- (((input_df$score - min(input_df$score)) / (max(input_df$score) - min(input_df$score))) * (4)) + 1

        if (is.null(cell_cols)) {
            cell_cols <- setNames(paletteMartin(n = length(unique(c(input_df$source, input_df$target)))), unique(c(input_df$source, input_df$target)))
        }

        segments <- unique(c(paste0(input_df$source, "|", input_df$ligand), paste0(input_df$target, "|", input_df$receptor)))
        grp <- str_extract(segments, "[^|]+")
        names(grp) <- segments

        gene_df <- as.data.frame(exp_df %>% mutate(cell_gene = paste0(cell_type, "|", gene)) %>% filter(cell_gene %in% segments))
        rownames(gene_df) <- gene_df$cell_gene

        brks <- scales::pretty_breaks(n = 5)(c(floor(min(gene_df$mean_exp)), ceiling(max(gene_df$mean_exp))))
        gene_col_fun <- colorRamp2(brks, RColorBrewer::brewer.pal(length(brks), palette))

        inner.cols <- setNames(gene_col_fun(gene_df[segments, "mean_exp"]), segments)
        lgd1 <- Legend(
            labels = unique(c(input_df$source, input_df$target)),
            title = "Cell type",
            type = "points",
            title_gp = gpar(fontsize = 14 * cex),
            labels_gp = gpar(fontsize = 12 * cex),
            legend_gp = gpar(col = "transparent"),
            background = cell_cols[unique(c(input_df$source, input_df$target))],
            direction = "horizontal"
        )

        lgd2 <- Legend(
            title_gp = gpar(fontsize = 14 * cex),
            labels_gp = gpar(fontsize = 12 * cex),
            direction = "horizontal", at = brks,
            col_fun = gene_col_fun, title = "Mean exp."
        )
        circlize_plot <- function() {
            par(cex = cex)
            chordDiagram(
                input_df %>%
                    select(source_lig, target_rec, score),
                directional = 1, group = grp, link.sort = FALSE, diffHeight = 0.005, scale = scale,
                direction.type = c("arrows"), link.arr.type = "triangle", annotationTrack = c(),
                preAllocateTracks = list(list(track.height = 0.175), list(track.height = 0.05), list(track.height = 0.045)),
                big.gap = 3, transparency = 1, link.arr.lwd = arr_wd, link.arr.col = "black", link.arr.length = 0.4, link.arr.width = 0.35, ...
            )
            circos.track(track.index = 1, panel.fun = function(x, y) {
                circos.text(CELL_META$xcenter, CELL_META$ylim[1], str_extract(CELL_META$sector.index, "[^|]+$"),
                    facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1.3
                )
            }, bg.border = NA)
            for (l in unique(str_extract(segments, "[^|]+"))) {
                highlight.sector(segments[str_detect(segments, paste0("^", str_escape(l)))], track.index = 2, col = cell_cols[l])
            }
            circos.track(track.index = 3, panel.fun = function(x, y) {
                circos.rect(CELL_META$xlim[1], CELL_META$ylim[1], CELL_META$xlim[2], CELL_META$ylim[2],
                    sector.index = CELL_META$sector.index, col = inner.cols[CELL_META$sector.index]
                )
            }, bg.border = NA)
            if (show_legend == TRUE) {
                draw(packLegend(lgd1, lgd2, direction = "vertical"), just = c("left", "bottom"), x = unit(4.75, "mm"), y = unit(4.75, "mm"))
            }
            circos.clear()
        }
    }
    circlize_plot()
}
