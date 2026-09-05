############## downstream ##############
############## DMEs ##############
if (length(dmes) > 0) {
    dmes <- lapply(dmes, function(x) list_update(dmes_defaults, x))
    dmes_dir <- file.path(tables_dir, "dmes")
    dir.create(dmes_dir, recursive = TRUE, showWarnings = FALSE)
    for (name in names(dmes)) {
        case <- dmes[[name]]
        tryCatch({
            log$info(glue("Running DME analysis `{name}` ..."))
            mode <- case$mode %||% "find_all"
            if (!mode %in% c("find", "find_all")) {
                stop(glue("Unknown DME mode `{mode}` for `{name}`!"))
            }
            if (mode == "find" && (is.null(case$barcodes1) || is.null(case$barcodes2))) {
                stop("`barcodes1` and `barcodes2` are required for `mode = \"find\"`!")
            }
            info <- case_info(paste0("DMEs::", name), dmes_dir, is_dir = FALSE, create = TRUE)
            lollipop <- case$lollipop %||% TRUE
            volcano <- case$volcano %||% FALSE
            devpars <- case$devpars %||% list()
            more_formats <- case$more_formats %||% list()
            descr <- case$descr %||% "The results of the differential module expression analysis."
            fn_args <- case[setdiff(names(case), c("mode", "lollipop", "volcano", "devpars", "more_formats", "descr", "pvalue", "plot_labels", "show_cutoff"))]
            if (mode == "find") {
                dmes_res <- do_call("FindDMEs", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), fn_args))
            } else {
                dmes_res <- do_call("FindAllDMEs", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), fn_args))
            }
            write_table(dmes_res, paste0(info$prefix, ".tsv"))
            report_items <- list(list(kind = "descr", content = descr))
            # PlotDMEsLollipop's group.by must be a column in the DMEs table
            # (for FindAllDMEs it is `group`); map it there instead of dropping
            # it, otherwise PlotDMEsLollipop merges all comparisons and
            # PlotLollipop fails on duplicated module levels
            plot_case <- case
            if (mode == "find") {
                # two-group DMEs (FindDMEs) have no comparison column;
                # hdWGCNA's plotters treat a missing group.by as "all in
                # one group" (the tutorial pattern), while an explicit
                # group.by that is absent from the DMEs data errors out
                plot_case$group.by <- NULL
            } else if (!is.null(plot_case$group.by) && !(plot_case$group.by %in% names(dmes_res)) && "group" %in% names(dmes_res)) {
                plot_case$group.by <- "group"
            }
            plot_case$pvalue <- plot_case$pvalue %||% "p_val_adj"
            # the plotters return a list of per-group plots when there are
            # multiple groups (FindAllDMEs); wrap them into one tall stack
            # and size the device to the number of groups
            if (lollipop) {
                p <- do_call("PlotDMEsLollipop", c(
                    list(seurat_obj = srtobj, DMEs = dmes_res, wgcna_name = wgcna_name),
                    hd_args(plot_case, "PlotDMEsLollipop")
                ))
                if (is.list(p) && !inherits(p, "ggplot")) {
                    save_plot_list(p, paste0(info$prefix, ".lollipop"), auto_devpars(devpars, length(p), ncol = 1), ncol = 1, formats = "png")
                } else {
                    save_plot(p, paste0(info$prefix, ".lollipop"), devpars, formats = "png")
                }
                report_items <- c(report_items, list(reporter$image(paste0(info$prefix, ".lollipop"), more_formats, FALSE, kind = "image")))
            }
            if (volcano) {
                p <- do_call("PlotDMEsVolcano", c(
                    list(seurat_obj = srtobj, DMEs = dmes_res, wgcna_name = wgcna_name),
                    hd_args(plot_case, "PlotDMEsVolcano")
                ))
                if (is.list(p) && !inherits(p, "ggplot")) {
                    save_plot_list(p, paste0(info$prefix, ".volcano"), auto_devpars(devpars, length(p), ncol = 1), ncol = 1, formats = "png")
                } else {
                    save_plot(p, paste0(info$prefix, ".volcano"), devpars, formats = "png")
                }
                report_items <- c(report_items, list(reporter$image(paste0(info$prefix, ".volcano"), more_formats, FALSE, kind = "image")))
            }
            reporter$add2(
                list(name = "Plot", contents = report_items),
                list(name = "Table", contents = list(
                    list(kind = "descr", content = "The table of the DME analysis results."),
                    list(kind = "table", src = paste0(info$prefix, ".tsv"), data = list(nrows = 100))
                )),
                hs = hs_section("DMEs", name),
                ui = "tabs"
            )
        }, error = function(e) {
            log$warn(glue("Failed to run DME analysis `{name}`: {conditionMessage(e)}"))
        })
    }
}

############## Module-trait correlation ##############
if (length(module_trait_corr) > 0) {
    module_trait_corr <- lapply(module_trait_corr, function(x) list_update(module_trait_corr_defaults, x))
    mtc_dir <- file.path(tables_dir, "module_trait_corr")
    dir.create(mtc_dir, recursive = TRUE, showWarnings = FALSE)
    for (name in names(module_trait_corr)) {
        case <- module_trait_corr[[name]]
        tryCatch({
            log$info(glue("Running module-trait correlation `{name}` ..."))
            if (is.null(case$traits)) stop("`traits` is required for module-trait correlation!")
            info <- case_info(paste0("Module-Trait Correlation::", name), mtc_dir, is_dir = FALSE, create = TRUE)
            devpars <- case$devpars %||% list()
            more_formats <- case$more_formats %||% list()
            descr <- case$descr %||% "The correlation between the module eigengenes and the traits."
            mtc_args <- case[setdiff(names(case), c("devpars", "more_formats", "descr", "plot_args"))]
            srtobj <<- do_call("ModuleTraitCorrelation", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), mtc_args))
            mtc <- GetModuleTraitCorrelation(srtobj, wgcna_name = wgcna_name)
            for (k in names(mtc)) {
                df_list <- mtc[[k]]
                # the per-group matrices must become data.frames first:
                # `$<-` on a matrix flattens it into a list
                df_all <- do.call(rbind, Map(function(df, g) {
                    df <- as.data.frame(df)
                    df$group <- g
                    df
                }, df_list, names(df_list)))
                write_table(df_all, paste0(info$prefix, ".", k, ".tsv"), row.names = TRUE)
            }
            ps <- do_call("PlotModuleTraitCorrelation", c(
                list(seurat_obj = srtobj, wgcna_name = wgcna_name),
                list_update(hd_args(case, "PlotModuleTraitCorrelation"), case$plot_args %||% list())
            ))
            save_plot_list(ps, info$prefix, auto_devpars(devpars, length(ps), ncol = 1), ncol = 1, formats = c("png", more_formats))
            reporter$add2(
                list(name = "Plot", contents = list(
                    list(kind = "descr", content = descr),
                    reporter$image(info$prefix, more_formats, FALSE, kind = "image")
                )),
                list(name = "Correlation", contents = list(
                    list(kind = "descr", content = "The correlation between the module eigengenes and the traits, per group."),
                    list(kind = "table", src = paste0(info$prefix, ".cor.tsv"), data = list(nrows = 100))
                )),
                list(name = "P-value", contents = list(
                    list(kind = "table", src = paste0(info$prefix, ".pval.tsv"), data = list(nrows = 100))
                )),
                list(name = "FDR", contents = list(
                    list(kind = "table", src = paste0(info$prefix, ".fdr.tsv"), data = list(nrows = 100))
                )),
                hs = hs_section("Module-Trait Correlation", name),
                ui = "tabs"
            )
        }, error = function(e) {
            log$warn(glue("Failed to run module-trait correlation `{name}`: {conditionMessage(e)}"))
        })
    }
}

############## Enrichment ##############
if (length(enrich) > 0) {
    enrich <- lapply(enrich, function(x) list_update(enrich_defaults, x))
    enrich_dir <- file.path(tables_dir, "enrich")
    enrich_plots_dir <- file.path(plots_dir, "Enrich")
    dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(enrich_plots_dir, recursive = TRUE, showWarnings = FALSE)
    mod_table <- do_call("GetModules", list(seurat_obj = srtobj, wgcna_name = wgcna_name))
    # `module` is a factor — compare labels, not level codes (see HdWGCNA-plots.R)
    mod_table <- mod_table[as.character(mod_table$module) != "grey", , drop = FALSE]
    mods <- split(mod_table$gene_name, mod_table$module)
    for (name in names(enrich)) {
        case <- enrich[[name]]
        tryCatch({
            log$info(glue("Running enrichment analysis `{name}` ..."))
            if (is.null(case$dbs)) stop("`dbs` is required for enrichment analysis!")
            info <- case_info(paste0("Enrichment::", name), enrich_dir, is_dir = FALSE, create = TRUE)
            devpars <- case$devpars %||% list()
            more_formats <- case$more_formats %||% list()
            descr <- case$descr %||% glue("The functional enrichment analysis of the module genes (dbs: {paste(case$dbs, collapse = ', ')}).")
            enrich_plots <- case$enrich_plots %||% list()
            dbs <- case$dbs
            for (mod in names(mods)) {
                genes <- mods[[mod]]
                if (length(genes) < 5) {
                    log$warn(glue("Module `{mod}` has fewer than 5 genes, skipped for enrichment in `{name}`."))
                    next
                }
                enrich_df <- RunEnrichment(genes, dbs = dbs, style = "enrichr")
                if (is.null(enrich_df) || nrow(enrich_df) == 0) {
                    log$warn(glue("No enrichment results for module `{mod}` in `{name}`."))
                    next
                }
                write_table(enrich_df, paste0(info$prefix, ".", mod, ".tsv"))
                mod_items <- list(list(kind = "descr", content = glue("Enrichment analysis for module `{mod}` ({length(genes)} genes).")))
                for (db in unique(enrich_df$Database)) {
                    for (plotname in names(enrich_plots)) {
                        plotargs <- extract_vars(enrich_plots[[plotname]], "descr", allow_nonexisting = TRUE)
                        if (!is.null(plotargs$db) && !identical(plotargs$db, db)) next
                        plotargs$data <- enrich_df[enrich_df$Database == db, , drop = FALSE]
                        p <- do_call(VizEnrichment, plotargs)
                        outprefix <- file.path(enrich_plots_dir, paste0(name, ".", mod, ".", slugify(db), ".", slugify(plotname)))
                        save_plot(p, outprefix, plotargs$devpars, formats = "png")
                        mod_items <- c(mod_items, list(reporter$image(outprefix, NULL, FALSE, kind = "image")))
                    }
                }
                mod_items <- c(mod_items, list(list(kind = "table", src = paste0(info$prefix, ".", mod, ".tsv"), data = list(nrows = 100))))
                reporter$add2(
                    list(name = mod, contents = mod_items),
                    hs = hs_section("Enrichment", name),
                    ui = "tabs"
                )
            }
        }, error = function(e) {
            log$warn(glue("Failed to run enrichment analysis `{name}`: {conditionMessage(e)}"))
        })
    }
}

############## GSEA ##############
if (length(gsea) > 0) {
    gsea <- lapply(gsea, function(x) list_update(gsea_defaults, x))
    gsea_dir <- file.path(tables_dir, "gsea")
    gsea_plots_dir <- file.path(plots_dir, "GSEA")
    dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(gsea_plots_dir, recursive = TRUE, showWarnings = FALSE)
    mod_table <- do_call("GetModules", list(seurat_obj = srtobj, wgcna_name = wgcna_name))
    kme_cols <- grep("^kME_", colnames(mod_table), value = TRUE)
    if (length(kme_cols) == 0) {
        log$warn("No kME columns in the modules table, skipping GSEA ...")
    } else {
        for (name in names(gsea)) {
            case <- gsea[[name]]
            tryCatch({
                log$info(glue("Running GSEA `{name}` ..."))
                genesets <- case$genesets
                if (is.null(genesets)) stop("`genesets` is required for GSEA!")
                if (is.list(genesets)) {
                    # arguments for `msigdbr::msigdbr()`, e.g.
                    # {species: "human", collection: "C5", subcollection: "BP"}
                    gs_df <- do.call(msigdbr::msigdbr, genesets)
                    genesets <- split(gs_df$gene_symbol, gs_df$gs_name)
                    genesets <- lapply(genesets, unique)
                } else if (file.exists(genesets)) {
                    genesets <- ParseGMT(genesets)
                } else {
                    stop(glue("`genesets` must be a GMT file or a dict of arguments for `msigdbr::msigdbr()`!"))
                }
                metric <- case$metric %||% "p.adjust"
                cutoff <- case$cutoff %||% 0.05
                top_term <- case$top_term %||% 10
                descr <- case$descr %||% "The gene set enrichment analysis of the module gene ranks (kME) against the gene sets."
                info <- case_info(paste0("GSEA::", name), gsea_dir, is_dir = FALSE, create = TRUE)
                more_formats <- case$more_formats %||% list()
                devpars <- case$devpars %||% list()
                gsea_plots <- case$gsea_plots %||% list()
                # fgsea output columns are pval/padj; plotthis uses pvalue/p.adjust
                metric_col <- if (metric == "p.adjust") "padj" else if (metric == "pvalue") "pval" else metric
                rest <- case[setdiff(names(case), c("descr", "devpars", "more_formats", "metric", "cutoff", "top_term", "gsea_plots", "genesets"))]
                for (mod in setdiff(unique(as.character(mod_table$module)), "grey")) {
                    kme_col <- paste0("kME_", mod)
                    if (!kme_col %in% colnames(mod_table)) {
                        log$warn(glue("No kME column `{kme_col}` for module `{mod}`, skipped in GSEA `{name}`."))
                        next
                    }
                    ranks <- setNames(mod_table[[kme_col]], mod_table$gene_name)
                    ranks <- ranks[!is.na(ranks)]
                    if (length(ranks) < 5) {
                        log$warn(glue("Module `{mod}` has fewer than 5 genes with kME values, skipped in GSEA `{name}`."))
                        next
                    }
                    log$info(glue("Running GSEA for module `{mod}` ..."))
                    result <- do_call("RunGSEA", c(list(ranks = ranks, genesets = genesets), rest))
                    write_table(result, paste0(info$prefix, ".", mod, ".tsv"))
                    mod_items <- list(list(kind = "descr", content = glue("GSEA for module `{mod}` ({length(ranks)} genes ranked by kME).")))
                    # summary dot plot (VizGSEA's first argument is
                    # `gsea_results`, not `data`)
                    summary_args <- list(gsea_results = result, plot_type = "summary", top_term = top_term, signif_by = metric, signif_cutoff = cutoff)
                    p <- do_call("VizGSEA", summary_args)
                    outprefix <- file.path(gsea_plots_dir, paste0(name, ".", mod, ".summary"))
                    save_plot(p, outprefix, devpars, formats = c("png", more_formats))
                    mod_items <- c(mod_items, list(reporter$image(outprefix, more_formats, FALSE, kind = "image")))
                    # running-score plots
                    for (plotname in names(gsea_plots)) {
                        plotargs <- gsea_plots[[plotname]]
                        # an empty case dict (`{}`) has no names; extract_vars
                        # requires a named list
                        if (is.null(names(plotargs))) names(plotargs) <- character(0)
                        plot_devpars <- plotargs$devpars %||% devpars
                        plotargs <- extract_vars(plotargs, "descr", "devpars", allow_nonexisting = TRUE)
                        gs <- plotargs$gs
                        if (is.null(gs)) {
                            # fgsea may leave pval/padj as NA for some gene
                            # sets; NA comparisons would leak NA into `gs`
                            gs <- head(result$pathway[!is.na(result[[metric_col]]) & result[[metric_col]] < cutoff], top_term)
                        }
                        if (length(gs) == 0) {
                            log$warn(glue("No significant terms for module `{mod}` in GSEA plot `{plotname}`."))
                            next
                        }
                        plotargs$gsea_results <- result
                        plotargs$plot_type <- "gsea"
                        plotargs$gs <- gs
                        plotargs$gene_sets <- attr(result, "gene_sets")[gs]
                        p <- do_call("VizGSEA", plotargs)
                        outprefix <- file.path(gsea_plots_dir, paste0(name, ".", mod, ".", slugify(plotname)))
                        save_plot(p, outprefix, plot_devpars, formats = c("png", more_formats))
                        mod_items <- c(mod_items, list(reporter$image(outprefix, more_formats, FALSE, kind = "image")))
                    }
                    mod_items <- c(mod_items, list(list(kind = "table", src = paste0(info$prefix, ".", mod, ".tsv"), data = list(nrows = 100))))
                    reporter$add2(
                        list(name = mod, contents = mod_items),
                        hs = hs_section("GSEA", name),
                        ui = "tabs"
                    )
                }
            }, error = function(e) {
                log$warn(glue("Failed to run GSEA `{name}`: {conditionMessage(e)}"))
            })
        }
    }
}
