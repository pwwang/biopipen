############## plots ##############
plot_sections <- c(
    soft_powers = "Soft Powers",
    dendrogram = "Dendrogram",
    kmes = "Module KMEs",
    module_umap = "Module UMAP",
    module_network = "Module Networks",
    hub_gene_network = "Hub Gene Networks",
    module_feature = "Module Features",
    module_radar = "Module Radars",
    module_correlogram = "Module Correlograms",
    consensus_compare = "Dendrogram"
)

default_descr <- c(
    soft_powers = "The scale-free topology model fit (SFT) and the mean connectivity for the tested soft powers. The soft power is chosen as the lowest power with an SFT R-squared of at least 0.8.",
    dendrogram = "The gene dendrogram with the detected modules colored.",
    kmes = "The kME (module membership) of the genes in each module.",
    module_umap = "The UMAP plot of the modules, colored by the module assignments.",
    module_network = "The module network plots showing the connectivity between the genes in each module.",
    hub_gene_network = "The hub gene network plot showing the connectivity between the hub genes of the modules.",
    module_feature = "The module feature plots showing the expression of the module eigengenes or scores.",
    module_radar = "The module radar plots showing the expression of the module eigengenes or scores across the groups.",
    module_correlogram = "The module correlogram showing the correlation between the module eigengenes.",
    consensus_compare = "The consensus dendrogram with the module colors of the consensus and the standard (non-consensus) workflows compared, as in the consensus WGCNA tutorial."
)

plots_dir <- file.path(outdir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

plots <- lapply(plots, function(x) list_update(plots_defaults, x))

# number of non-grey modules, for sizing the per-module plots (feature
# plots and radars wrap the per-module panels internally with ncol = 4)
# NOTE: `module` is a factor, so compare on character labels — `!=` on a
# factor compares level codes and counts every gene (grey is the last
# level) when "grey" is not the first level
n_modules <- tryCatch({
    mod_tbl <- do_call("GetModules", list(seurat_obj = srtobj, wgcna_name = wgcna_name))
    sum(as.character(mod_tbl$module) != "grey")
}, error = function(e) 1)

# Some plotters (dendrogram, module UMAP, hub gene network, correlogram) draw
# on the current device via base graphics / corrplot / igraph and return NULL
# or an invisible list. save_plot would print that and produce a blank image,
# so open the device first, let the plotter draw, then close it.
save_plot_sidefx <- function(expr, prefix, devpars, formats = "png") {
    devpars <- list_update(list(res = 100, width = 8, height = 6), devpars)
    for (format in formats) {
        file <- paste0(prefix, ".", format)
        if (format == "png") {
            grDevices::png(file, width = devpars$width, height = devpars$height,
                units = "in", res = devpars$res)
        } else if (format == "pdf") {
            grDevices::pdf(file, width = devpars$width, height = devpars$height)
        } else {
            stop(glue("Unsupported format `{format}` for side-effect plot!"))
        }
        tryCatch(force(expr), finally = grDevices::dev.off())
    }
}

# ModuleUMAPPlot and HubGeneNetworkPlot color the edge data frame by
# iterating over ~1M rows with per-row `[.data.frame` extraction. That loop
# segfaults nondeterministically inside R's own subsetting code under the
# huge captured frame (the whole Seurat object). Replace the closure with
# the equivalent vectorized match()/ifelse computation (verified to give
# identical colors on real data). Idempotent: skips if already patched.
vectorize_edge_color <- function(fname, repl) {
    ns <- asNamespace("hdWGCNA")
    f <- get(fname, ns)
    b <- body(f)
    patched <- FALSE
    for (i in seq_along(b)) {
        if (grepl("future_sapply", paste(deparse(b[[i]]), collapse = " "))) {
            b[[i]] <- repl
            patched <- TRUE
            break
        }
    }
    if (!patched) return(invisible(NULL))
    body(f) <- b
    if (bindingIsLocked(fname, ns)) {
        unlockBinding(fname, ns)
    }
    assign(fname, f, envir = ns)
    # symbol lookup (do_call by name) resolves via the attached package env
    penv <- as.environment("package:hdWGCNA")
    if (bindingIsLocked(fname, penv)) {
        unlockBinding(fname, penv)
    }
    assign(fname, f, envir = penv)
    log$info(glue("Vectorized the edge coloring in hdWGCNA::{fname}"))
}
vectorize_edge_color("ModuleUMAPPlot", quote(edge_df$color <- ifelse(
    selected_modules$color[match(edge_df$Var1, selected_modules$gene_name)] ==
        selected_modules$color[match(edge_df$Var2, selected_modules$gene_name)],
    selected_modules$color[match(edge_df$Var1, selected_modules$gene_name)],
    "grey90"
)))
vectorize_edge_color("HubGeneNetworkPlot", quote(edge_df$color <- ifelse(
    modules$color[match(edge_df$Var1, modules$gene_name)] ==
        modules$color[match(edge_df$Var2, modules$gene_name)],
    modules$color[match(edge_df$Var1, modules$gene_name)],
    "grey90"
)))

# module umap plots first, since they modify the object
for (name in names(plots)) {
    case <- plots[[name]]
    if (is.null(case$kind) || !identical(case$kind, "module_umap")) next
    tryCatch({
        log$info(glue("Plotting module UMAP for `{name}` ..."))
        srtobj <<- do_call("RunModuleUMAP", c(
            list(seurat_obj = srtobj, wgcna_name = wgcna_name),
            case[setdiff(names(case), c("kind", "devpars", "more_formats", "descr", "umap_plot_args"))]
        ))
        umap_args <- list_update(list(), case$umap_plot_args %||% list())
        info <- case_info(paste0("Module UMAP::", name), plots_dir, is_dir = FALSE, create = TRUE)
        save_plot_sidefx(
            do_call("ModuleUMAPPlot", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), umap_args)),
            info$prefix, case$devpars, formats = c("png", case$more_formats)
        )
        reporter$add2(
            list(name = "Plot", contents = list(
                list(kind = "descr", content = case$descr %||% default_descr[["module_umap"]]),
                reporter$image(info$prefix, case$more_formats, FALSE, kind = "image")
            )),
            list(name = "Hub Genes", contents = list(
                list(kind = "descr", content = "The hub genes of the modules, which are labeled on the UMAP plot."),
                list(kind = "table", src = file.path(tables_dir, "hub_genes.tsv"), data = list(nrows = 100))
            )),
            hs = hs_section("Module UMAP", name),
            ui = "tabs"
        )
    }, error = function(e) {
        log$warn(glue("Failed to plot module UMAP for `{name}`: {conditionMessage(e)}"))
    })
}

for (name in names(plots)) {
    case <- plots[[name]]
    kind <- case$kind
    if (is.null(kind) || identical(kind, "module_umap")) next
    tryCatch({
        log$info(glue("Plotting `{kind}` for `{name}` ..."))
        section <- plot_sections[[kind]] %||% "Plots"
        info <- case_info(paste0(section, "::", name), plots_dir, is_dir = FALSE, create = TRUE)
        plot_args <- case[setdiff(names(case), c("kind", "devpars", "more_formats", "descr"))]
        switch(
            kind,
            soft_powers = {
                ps <- do_call("PlotSoftPowers", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), plot_args))
                # non-consensus: unnamed list of 1-4 ggplots (SFT fit, mean, median, max
                # connectivity); consensus: list of per-set lists of ggplots (group-major,
                # as in the consensus tutorial: plot_list[[i]][[1]] per set). Transpose
                # the nested list so each output shows one plot type across all sets
                # side-by-side; save_plot on a plain list draws every element onto the
                # same device and leaves only the last one visible, hence the wrap.
                suffixes <- c("sft", "fit", "median.k", "max.k")
                ps_types <- if (length(ps) > 0 && is.list(ps[[1]]) && !inherits(ps[[1]], "ggplot")) {
                    lapply(seq_along(suffixes), function(j) lapply(ps, function(g) g[[j]]))
                } else {
                    ps
                }
                for (i in seq_along(ps_types)) {
                    this <- ps_types[[i]]
                    if (is.list(this) && !inherits(this, "ggplot")) {
                        this <- patchwork::wrap_plots(this, ncol = length(this))
                    }
                    save_plot(this, paste0(info$prefix, ".", suffixes[[i]]), auto_devpars(case$devpars, length(ps_types[[i]])), formats = c("png", case$more_formats))
                }
                reporter$add2(
                    list(name = "Plot", contents = list(
                        list(kind = "descr", content = case$descr %||% default_descr[["soft_powers"]]),
                        reporter$image(paste0(info$prefix, ".sft"), case$more_formats, FALSE, kind = "image"),
                        reporter$image(paste0(info$prefix, ".fit"), case$more_formats, FALSE, kind = "image"),
                        reporter$image(paste0(info$prefix, ".median.k"), case$more_formats, FALSE, kind = "image"),
                        reporter$image(paste0(info$prefix, ".max.k"), case$more_formats, FALSE, kind = "image")
                    )),
                    hs = hs_section("Introduction", name),
                    ui = "tabs"
                )
            },
            dendrogram = {
                save_plot_sidefx(
                    do_call("PlotDendrogram", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), plot_args)),
                    info$prefix, case$devpars, formats = c("png", case$more_formats)
                )
                reporter$add2(
                    list(name = "Plot", contents = list(
                        list(kind = "descr", content = case$descr %||% default_descr[["dendrogram"]]),
                        reporter$image(info$prefix, case$more_formats, FALSE, kind = "image")
                    )),
                    hs = hs_section("Introduction", name),
                    ui = "tabs"
                )
            },
            consensus_compare = {
                # the consensus tutorial's comparison figure: re-run the
                # standard (non-consensus) workflow on the same data under a
                # different wgcna_name, then plot both module-color tracks on
                # the consensus dendrogram
                if (!use_consensus) {
                    stop("`consensus_compare` plot kind requires `use_consensus`!")
                }
                standard_name <- plot_args$standard_name
                if (is.null(standard_name)) {
                    stop("`consensus_compare` plot kind requires `standard_name`!")
                }
                soft_power <- plot_args$soft_power
                main <- plot_args$main %||% "Consensus vs standard dendrogram"
                test_soft_powers <- is.null(soft_power)
                plot_args$standard_name <- NULL
                plot_args$soft_power <- NULL
                plot_args$main <- NULL
                # `wgcna_name` is passed explicitly in the do_call below; adding
                # it to plot_args would duplicate the argument
                # run the standard workflow with the same setup and data
                standard_setup_args <- SetupForWGCNAArgs
                standard_setup_args$wgcna_name <- NULL
                srtobj <- do_call("SetupForWGCNA", c(list(seurat_obj = srtobj, wgcna_name = standard_name), standard_setup_args))
                if (use_pseudobulk) {
                    plot_args$mat <- se
                    if (is.null(plot_args$layer)) plot_args$layer <- "VST"
                    plot_args$use_metacells <- NULL
                    plot_args$group_name <- NULL
                    plot_args$group.by <- NULL
                }
                log$info(glue("Setting up the standard (non-consensus) expression data for `{standard_name}` ..."))
                srtobj <- do_call("SetDatExpr", c(list(seurat_obj = srtobj, wgcna_name = standard_name), plot_args))
                if (test_soft_powers) {
                    log$info(glue("Testing the soft powers for `{standard_name}` ..."))
                    srtobj <- do_call("TestSoftPowers", list(seurat_obj = srtobj, wgcna_name = standard_name))
                }
                log$info(glue("Constructing the standard network `{standard_name}` ..."))
                # a standard (non-consensus) run: no consensus flag and a
                # distinct TOM name so the consensus TOM is not overwritten
                standard_net_args <- ConstructNetworkArgs
                standard_net_args$consensus <- FALSE
                standard_net_args$tom_name <- NULL
                if (!is.null(soft_power)) standard_net_args$soft_power <- soft_power
                srtobj <- do_call("ConstructNetwork", c(list(seurat_obj = srtobj, wgcna_name = standard_name), standard_net_args))
                log$info(glue("Computing the module eigengenes and connectivity for `{standard_name}` ..."))
                srtobj <- do_call("ModuleEigengenes", list(seurat_obj = srtobj, wgcna_name = standard_name))
                srtobj <- do_call("ModuleConnectivity", list(seurat_obj = srtobj, wgcna_name = standard_name))
                # build the color comparison on the consensus dendrogram
                consensus_modules <- do_call("GetModules", list(seurat_obj = srtobj, wgcna_name = wgcna_name))
                standard_modules <- do_call("GetModules", list(seurat_obj = srtobj, wgcna_name = standard_name))
                net <- srtobj@misc[[wgcna_name]]$wgcna_net
                dendro <- net$dendrograms[[1]]
                consensus_colors <- consensus_modules$color
                names(consensus_colors) <- consensus_modules$gene_name
                standard_colors <- standard_modules$color
                names(standard_colors) <- standard_modules$gene_name
                color_df <- data.frame(
                    consensus = consensus_colors,
                    standard = standard_colors[consensus_modules$gene_name],
                    row.names = NULL
                )
                save_plot_sidefx(
                    WGCNA::plotDendroAndColors(
                        dendro, color_df,
                        groupLabels = colnames(color_df),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = main
                    ),
                    info$prefix, case$devpars, formats = c("png", case$more_formats)
                )
                # the standard workflow dendrogram (tutorial figure 3)
                save_plot_sidefx(
                    do_call("PlotDendrogram", list(seurat_obj = srtobj, wgcna_name = standard_name, main = main)),
                    paste0(info$prefix, ".standard"), case$devpars, formats = c("png", case$more_formats)
                )
                reporter$add2(
                    list(name = "Plot", contents = list(
                        list(kind = "descr", content = case$descr %||% default_descr[["consensus_compare"]]),
                        reporter$image(info$prefix, case$more_formats, FALSE, kind = "image"),
                        reporter$image(paste0(info$prefix, ".standard"), case$more_formats, FALSE, kind = "image")
                    )),
                    hs = hs_section(section, name),
                    ui = "tabs"
                )
            },
            kmes = {
                ps <- do_call("PlotKMEs", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), plot_args))
                save_plot_list(ps, info$prefix, auto_devpars(case$devpars, length(ps), ncol = 5), ncol = 5, formats = c("png", case$more_formats))
                reporter$add2(
                    list(name = "Plot", contents = list(
                        list(kind = "descr", content = case$descr %||% default_descr[["kmes"]]),
                        reporter$image(info$prefix, case$more_formats, FALSE, kind = "image")
                    )),
                    hs = hs_section(section, name),
                    ui = "tabs"
                )
            },
            module_network = {
                plot_args$outdir <- file.path(plots_dir, paste0("ModuleNetworks.", info$slug))
                dir.create(plot_args$outdir, recursive = TRUE, showWarnings = FALSE)
                do_call("ModuleNetworkPlot", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), plot_args))
                # ModuleNetworkPlot writes one PDF per module; merge them into
                # a single PDF so the report embeds one file instead of dozens
                net_pdfs <- list.files(plot_args$outdir, pattern = "\\.pdf$", full.names = TRUE)
                net_items <- list(list(kind = "descr", content = case$descr %||% default_descr[["module_network"]]))
                merged <- file.path(plot_args$outdir, "module_networks.pdf")
                if (length(net_pdfs) > 0 && merge_pdfs(net_pdfs, merged)) {
                    net_items <- c(net_items, list(list(kind = "pdf", src = merged)))
                } else {
                    net_items <- c(net_items, lapply(net_pdfs, function(f) list(kind = "pdf", src = f)))
                }
                reporter$add2(
                    list(name = "Plot", contents = net_items),
                    hs = hs_section(section, name),
                    ui = "tabs"
                )
            },
            hub_gene_network = {
                save_plot_sidefx(
                    do_call("HubGeneNetworkPlot", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), plot_args)),
                    info$prefix, case$devpars, formats = c("png", case$more_formats)
                )
                reporter$add2(
                    list(name = "Plot", contents = list(
                        list(kind = "descr", content = case$descr %||% default_descr[["hub_gene_network"]]),
                        reporter$image(info$prefix, case$more_formats, FALSE, kind = "image")
                    )),
                    hs = hs_section(section, name),
                    ui = "tabs"
                )
            },
            module_feature = {
                ps <- do_call("ModuleFeaturePlot", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), plot_args))
                save_plot_list(ps, info$prefix, auto_devpars(case$devpars, length(ps), ncol = 4), ncol = 4, formats = c("png", case$more_formats))
                reporter$add2(
                    list(name = "Plot", contents = list(
                        list(kind = "descr", content = case$descr %||% default_descr[["module_feature"]]),
                        reporter$image(info$prefix, case$more_formats, FALSE, kind = "image")
                    )),
                    hs = hs_section(section, name),
                    ui = "tabs"
                )
            },
            module_radar = {
                ps <- do_call("ModuleRadarPlot", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), plot_args))
                save_plot_list(ps, info$prefix, auto_devpars(case$devpars, length(ps), ncol = 4), ncol = 4, formats = c("png", case$more_formats))
                reporter$add2(
                    list(name = "Plot", contents = list(
                        list(kind = "descr", content = case$descr %||% default_descr[["module_radar"]]),
                        reporter$image(info$prefix, case$more_formats, FALSE, kind = "image")
                    )),
                    hs = hs_section(section, name),
                    ui = "tabs"
                )
            },
            module_correlogram = {
                save_plot_sidefx(
                    do_call("ModuleCorrelogram", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), plot_args)),
                    info$prefix, case$devpars, formats = c("png", case$more_formats)
                )
                reporter$add2(
                    list(name = "Plot", contents = list(
                        list(kind = "descr", content = case$descr %||% default_descr[["module_correlogram"]]),
                        reporter$image(info$prefix, case$more_formats, FALSE, kind = "image")
                    )),
                    hs = hs_section(section, name),
                    ui = "tabs"
                )
            },
            stop(glue("Unknown plot kind `{kind}` for `{name}`!"))
        )
    }, error = function(e) {
        log$warn(glue("Failed to plot `{name}`: {conditionMessage(e)}"))
    })
}
