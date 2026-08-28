############## motifs ##############
# Replicates hdWGCNA's MotifOverlapBarPlot (broken in 0.4.12, see below).
motif_overlap_bar_plot <- function(
    seurat_obj, wgcna_name, outdir,
    n_tfs = 10, plot_size = c(5, 6),
    motif_font = "helvetica_regular", module_names = NULL
) {
    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }
    modules <- GetModules(seurat_obj)
    mods <- levels(modules$module)
    mods <- mods[mods != "grey"]
    if (is.null(module_names)) {
        module_names <- mods
    }
    overlap_df <- GetMotifOverlap(seurat_obj, wgcna_name)
    motif_df <- GetMotifs(seurat_obj)
    pfm <- GetPFMList(seurat_obj)
    overlap_df$motif_ID <- motif_df$motif_ID[match(overlap_df$tf, motif_df$motif_name)]
    overlap_df <- overlap_df %>% subset(module %in% module_names)
    for (cur_mod in module_names) {
        log$info(glue("Plotting motif overlaps for module {cur_mod} ..."))
        plot_df <- overlap_df %>% subset(module == cur_mod) %>%
            top_n(n_tfs, wt = odds_ratio) %>% arrange(desc(odds_ratio))
        p1 <- plot_df %>% ggplot(aes(y = reorder(tf, odds_ratio), fill = odds_ratio, x = odds_ratio)) +
            geom_bar(stat = "identity", width = 0.7) + NoLegend() +
            scale_fill_gradient(high = unique(plot_df$color), low = "grey90") +
            ylab("") + theme(
                axis.line.y = element_blank(),
                axis.text.y = element_blank(),
                plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
            )
        plot_list <- list()
        for (i in seq_len(nrow(plot_df))) {
            cur_id <- plot_df[i, "motif_ID"]
            cur_name <- plot_df[i, "tf"]
            plot_list[[cur_id]] <- ggplot() +
                ggseqlogo::geom_logo(pfm[[cur_id]]@profileMatrix, font = motif_font) +
                ggseqlogo::theme_logo() +
                xlab("") + ylab(cur_name) + theme(
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    axis.title.y = element_text(angle = 0),
                    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
                )
        }
        patch1 <- wrap_plots(plot_list, ncol = 1)
        outplot <- (patch1 | p1) + plot_layout(ncol = 2, widths = c(1, 2)) +
            plot_annotation(
                title = paste0("Motif overlaps with ", cur_mod),
                theme = theme(plot.title = element_text(hjust = 0.5))
            )
        pdf(file.path(outdir, paste0(cur_mod, "_motif_overlaps.pdf")),
            width = plot_size[1], height = plot_size[2], useDingbats = FALSE)
        print(outplot)
        dev.off()
    }
}

if (length(motifs) > 0) {
    tryCatch({
        if (!motif_scan_done) {
            stop("`tf_network.motif_scan` must be run before the motif overlap analysis!")
        }
        log$info("Computing the motif overlaps ...")
        srtobj <- do_call("OverlapModulesMotifs", list(seurat_obj = srtobj, wgcna_name = wgcna_name))

        motif_overlap <- GetMotifOverlap(srtobj, wgcna_name = wgcna_name)
        motif_df <- GetMotifs(srtobj)
        # MotifOverlapBarPlot resolves the PFM with
        # `motif_df$motif_ID[match(tf, motif_df$motif_name)]`, but MotifScan
        # keys the overlap `tf` by motif ID (e.g. MA0854.1), so the match
        # always fails and the plot dies on `pfm[[NA]]` (as.matrix(NULL)).
        # Map the IDs back to their names, and drop the overlaps that still
        # cannot be resolved (their motif's gene is absent from the object).
        ix <- match(motif_overlap$tf, motif_df$motif_ID)
        tf_sym <- motif_df$motif_name[ix]
        tf_sym[is.na(tf_sym)] <- motif_overlap$tf[is.na(tf_sym)]
        motif_overlap$tf <- tf_sym
        motif_overlap <- motif_overlap[!is.na(match(motif_overlap$tf, motif_df$motif_name)), ]
        srtobj <- do_call("SetMotifOverlap", list(seurat_obj = srtobj, overlap_df = motif_overlap, wgcna_name = wgcna_name))

        if (!is.null(motifs$overlap_bar_plot)) {
            obp <- motifs$overlap_bar_plot
            obp_dir <- file.path(plots_dir, "MotifOverlap")
            dir.create(obp_dir, recursive = TRUE, showWarnings = FALSE)
            obp_args <- obp[setdiff(names(obp), c("descr", "devpars", "more_formats"))]
            obp_args$outdir <- obp_dir
            log$info(glue("Plotting the motif overlaps to {obp_dir} ..."))
            # hdWGCNA's MotifOverlapBarPlot is broken in 0.4.12: the `proxy`
            # import shadows base `as.matrix` inside the hdWGCNA namespace, so
            # `as.matrix(pfm[[..]])` fails on PFMatrix objects. Replicate it
            # here using the `profileMatrix` slot directly (no dispatch).
            do.call(motif_overlap_bar_plot, c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), obp_args))
            descr <- obp$descr %||% "The motif overlap bar plots for the modules."
            obp_pdfs <- list.files(obp_dir, pattern = "_motif_overlaps\\.pdf$", full.names = TRUE)
            # merge the per-module PDFs into a single file so the report
            # embeds one PDF instead of dozens (fall back to the individual
            # files when ghostscript is unavailable)
            obp_items <- list(list(kind = "descr", content = descr))
            obp_merged <- file.path(obp_dir, "motif_overlaps.pdf")
            if (length(obp_pdfs) > 0 && merge_pdfs(obp_pdfs, obp_merged)) {
                obp_items <- c(obp_items, list(list(kind = "pdf", src = obp_merged)))
            } else {
                obp_items <- c(obp_items, lapply(obp_pdfs, function(f) list(kind = "pdf", src = f)))
            }
            reporter$add2(
                list(name = "Plot", contents = obp_items),
                list(name = "Table", contents = list(
                    list(kind = "descr", content = "The motif overlap table."),
                    list(kind = "table", src = file.path(tables_dir, "motif_overlap.tsv"), data = list(nrows = 100))
                )),
                hs = c("Motifs", "Motif Overlap"),
                ui = "tabs"
            )
        }

        motif_overlap <- GetMotifOverlap(srtobj, wgcna_name = wgcna_name)
        write_table(motif_overlap, file.path(tables_dir, "motif_overlap.tsv"))
    }, error = function(e) {
        log$warn(glue("Failed to run the motif overlap analysis: {conditionMessage(e)}"))
    })
}
