############## TF network ##############
motif_scan_done <- FALSE
if (length(tf_network) > 0) {
    tryCatch({
        tf_plots_dir <- file.path(plots_dir, "TF")
        dir.create(tf_plots_dir, recursive = TRUE, showWarnings = FALSE)

        ############## MotifScan ##############
        if (!is.null(tf_network$motif_scan)) {
            mscan <- tf_network$motif_scan
            if (is.null(mscan$species_genome)) stop("`species_genome` is required for MotifScan!")
            if (is.null(mscan$ensdb_package)) stop("`ensdb_package` is required for MotifScan!")
            log$info("Scanning the module genes for TF motifs ...")
            suppressPackageStartupMessages(library(mscan$ensdb_package, character.only = TRUE))
            jaspar_package <- mscan$jaspar_package %||% "JASPAR2020"
            suppressPackageStartupMessages(library(jaspar_package, character.only = TRUE))
            opts <- list(
                collection = mscan$collection %||% "CORE",
                tax_group = mscan$tax_group %||% "vertebrates",
                all_versions = isTRUE(mscan$all_versions)
            )
            pfm <- TFBSTools::getMatrixSet(x = get(jaspar_package), opts = opts)
            mscan_args <- mscan[setdiff(names(mscan), c("species_genome", "ensdb_package", "jaspar_package", "collection", "tax_group", "all_versions"))]
            srtobj <- do_call("MotifScan", c(list(
                seurat_obj = srtobj,
                pfm = pfm,
                EnsDb = get(mscan$ensdb_package),
                species_genome = mscan$species_genome,
                wgcna_name = wgcna_name
            ), mscan_args))
            motif_scan_done <- TRUE
        }

        ############## ConstructTFNetwork ##############
        if (!is.null(tf_network$construct)) {
            if (!motif_scan_done) {
                stop("`tf_network.motif_scan` must be set to construct the TF network!")
            }
            log$info("Expanding the gene set with the TFs and re-setting the expression data ...")
            motif_df <- GetMotifs(srtobj)
            tf_genes <- unique(motif_df$gene_name)
            mod_tbl <- GetModules(srtobj, wgcna_name = wgcna_name)
            # `module` is a factor — compare labels, not level codes (see HdWGCNA-plots.R)
            nongrey_genes <- mod_tbl$gene_name[as.character(mod_tbl$module) != "grey"]
            genes_use <- unique(c(tf_genes, nongrey_genes))
            srtobj <- do_call("SetWGCNAGenes", list(
                seurat_obj = srtobj,
                gene_list = genes_use,
                wgcna_name = wgcna_name
            ))
            srtobj <- do_call("SetDatExpr", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), SetDatExprArgs))
            model_params <- tf_network$construct$model_params
            if (is.null(model_params)) stop("`model_params` is required to construct the TF network!")
            if (is.null(model_params$nthread)) model_params$nthread <- ncores
            construct_args <- tf_network$construct[setdiff(names(tf_network$construct), "model_params")]
            construct_args$model_params <- model_params
            log$info("Constructing the TF network (xgboost) ...")
            srtobj <- do_call("ConstructTFNetwork", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), construct_args))
        }

        ############## AssignTFRegulons ##############
        if (!is.null(tf_network$assign)) {
            if (!motif_scan_done) {
                stop("`tf_network.motif_scan` must be run before assigning the TF regulons!")
            }
            log$info("Assigning the TF regulons ...")
            srtobj <- do_call("AssignTFRegulons", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), tf_network$assign))
        }

        ############## RegulonScores ##############
        if (!is.null(tf_network$regulon_scores)) {
            # FindDifferentialRegulons needs the scores for both target
            # types; run twice (the tutorial pattern) when target_type is
            # not set, once otherwise
            rscase <- tf_network$regulon_scores
            target_types <- rscase$target_type %||% c("positive", "negative")
            rscase$target_type <- NULL
            for (target_type in target_types) {
                log$info(glue("Computing the {target_type} regulon scores ..."))
                srtobj <- do_call("RegulonScores", c(
                    list(seurat_obj = srtobj, wgcna_name = wgcna_name, target_type = target_type),
                    rscase
                ))
            }
        }

        ############## TFNetworkPlot ##############
        if (length(tf_network$plots %||% list()) > 0) {
            for (pname in names(tf_network$plots)) {
                tryCatch({
                    tcase <- tf_network$plots[[pname]]
                    if (is.null(tcase$selected_tfs)) stop("`selected_tfs` is required for TFNetworkPlot!")
                    tinfo <- case_info(paste0("TF Network::", pname), tf_plots_dir, is_dir = FALSE, create = TRUE)
                    devpars <- tcase$devpars %||% list()
                    more_formats <- tcase$more_formats %||% list()
                    descr <- tcase$descr %||% glue("The TF regulatory network plot for {paste(tcase$selected_tfs, collapse = ', ')}.")
                    tcase <- tcase[setdiff(names(tcase), c("devpars", "more_formats", "descr"))]
                    p <- do_call("TFNetworkPlot", c(
                        list(seurat_obj = srtobj),
                        hd_args(c(list(wgcna_name = wgcna_name), tcase), "TFNetworkPlot")
                    ))
                    save_plot(p, tinfo$prefix, devpars, formats = c("png", more_formats))
                    reporter$add2(
                        list(name = "Plot", contents = list(
                            list(kind = "descr", content = descr),
                            reporter$image(tinfo$prefix, more_formats, FALSE, kind = "image")
                        )),
                        hs = hs_section("TF Network", pname),
                        ui = "tabs"
                    )
                }, error = function(e) {
                    log$warn(glue("Failed to plot TF network `{pname}`: {conditionMessage(e)}"))
                })
            }
        }

        ############## RegulonBarPlot ##############
        if (length(tf_network$regulon_bar_plots %||% list()) > 0) {
            for (pname in names(tf_network$regulon_bar_plots)) {
                tryCatch({
                    tcase <- tf_network$regulon_bar_plots[[pname]]
                    if (is.null(tcase$selected_tf)) stop("`selected_tf` is required for RegulonBarPlot!")
                    tinfo <- case_info(paste0("TF Network::", pname), tf_plots_dir, is_dir = FALSE, create = TRUE)
                    devpars <- tcase$devpars %||% list()
                    more_formats <- tcase$more_formats %||% list()
                    descr <- tcase$descr %||% glue("The regulon bar plot for TF {tcase$selected_tf}.")
                    tcase <- tcase[setdiff(names(tcase), c("devpars", "more_formats", "descr"))]
                    p <- do_call("RegulonBarPlot", c(
                        list(seurat_obj = srtobj),
                        hd_args(c(list(wgcna_name = wgcna_name), tcase), "RegulonBarPlot")
                    ))
                    save_plot(p, tinfo$prefix, devpars, formats = c("png", more_formats))
                    reporter$add2(
                        list(name = "Plot", contents = list(
                            list(kind = "descr", content = descr),
                            reporter$image(tinfo$prefix, more_formats, FALSE, kind = "image")
                        )),
                        hs = hs_section("TF Network", pname),
                        ui = "tabs"
                    )
                }, error = function(e) {
                    log$warn(glue("Failed to plot regulon bar plot `{pname}`: {conditionMessage(e)}"))
                })
            }
        }

        ############## FindDifferentialRegulons ##############
        # hdWGCNA 0.4.12 bug: the TF-expression DEG part of the function
        # references the bare variables `group1`/`group2` (line: degs <-
        # FindMarkers(exp_assay, cells.1 = group1, cells.2 = group2, ...)),
        # which are never defined in the function — the tutorial only works
        # because those globals happen to exist in the user's session. Patch
        # the namespace to use `barcodes1`/`barcodes2` (idempotent).
        ns <- asNamespace("hdWGCNA")
        f <- get("FindDifferentialRegulons", ns)
        b <- body(f)
        walk_dreg <- function(x) {
            if (is.symbol(x)) {
                if (identical(x, as.name("group1"))) return(as.name("barcodes1"))
                if (identical(x, as.name("group2"))) return(as.name("barcodes2"))
                return(x)
            }
            if (is.call(x) || is.pairlist(x)) {
                x[] <- lapply(x, walk_dreg)
            }
            x
        }
        b <- walk_dreg(b)
        body(f) <- b
        if (bindingIsLocked("FindDifferentialRegulons", ns)) {
            unlockBinding("FindDifferentialRegulons", ns)
        }
        assign("FindDifferentialRegulons", f, envir = ns)
        # do_call by name resolves via the attached package env, which still
        # points at the old function object — reassign there too
        if ("package:hdWGCNA" %in% search()) {
            penv <- as.environment("package:hdWGCNA")
            if (bindingIsLocked("FindDifferentialRegulons", penv)) {
                unlockBinding("FindDifferentialRegulons", penv)
            }
            assign("FindDifferentialRegulons", f, envir = penv)
        }
        log$info("Patched hdWGCNA::FindDifferentialRegulons (group1/group2 -> barcodes1/barcodes2)")
        dreg_dir <- file.path(tables_dir, "tf_network")
        dir.create(dreg_dir, recursive = TRUE, showWarnings = FALSE)
        if (length(tf_network$differential_regulons %||% list()) > 0) {
            for (dname in names(tf_network$differential_regulons)) {
                tryCatch({
                    dcase <- tf_network$differential_regulons[[dname]]
                    if (is.null(dcase$barcodes1) || is.null(dcase$barcodes2)) {
                        stop("`barcodes1` and `barcodes2` are required for FindDifferentialRegulons!")
                    }
                    dinfo <- case_info(paste0("TF Network::", dname), dreg_dir, is_dir = FALSE, create = TRUE)
                    devpars <- dcase$devpars %||% list()
                    more_formats <- dcase$more_formats %||% list()
                    descr <- dcase$descr %||% "The differential regulon analysis between the two groups of cells."
                    dcase <- dcase[setdiff(names(dcase), c("devpars", "more_formats", "descr"))]
                    dregs <- do_call("FindDifferentialRegulons", c(
                        list(seurat_obj = srtobj),
                        hd_args(c(list(wgcna_name = wgcna_name), dcase), "FindDifferentialRegulons")
                    ))
                    write_table(dregs, paste0(dinfo$prefix, ".tsv"))
                    p <- do_call("PlotDifferentialRegulons", c(
                        list(seurat_obj = srtobj, dregs = dregs),
                        hd_args(c(list(wgcna_name = wgcna_name), dcase), "PlotDifferentialRegulons")
                    ))
                    save_plot(p, dinfo$prefix, devpars, formats = c("png", more_formats))
                    reporter$add2(
                        list(name = "Plot", contents = list(
                            list(kind = "descr", content = descr),
                            reporter$image(dinfo$prefix, more_formats, FALSE, kind = "image")
                        )),
                        list(name = "Table", contents = list(
                            list(kind = "descr", content = "The results of the differential regulon analysis."),
                            list(kind = "table", src = paste0(dinfo$prefix, ".tsv"), data = list(nrows = 100))
                        )),
                        hs = hs_section("TF Network", dname),
                        ui = "tabs"
                    )
                }, error = function(e) {
                    log$warn(glue("Failed to run differential regulon analysis `{dname}`: {conditionMessage(e)}"))
                })
            }
        }

        ############## Regulon enrichment ##############
        # hdWGCNA's RunEnrichrRegulons uses the enrichR package; replicate its
        # flow with biopipen.utils::RunEnrichment (enrichit) and store the
        # results in the object so GetEnrichrRegulonTable still works below
        if (!is.null(tf_network$enrich_regulons)) {
            log$info("Running the regulon enrichment analysis ...")
            ereg <- tf_network$enrich_regulons
            dbs <- ereg$dbs %||% c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021")
            depth <- ereg$depth %||% 1
            use_regulons <- ereg$use_regulons %||% TRUE
            min_genes <- ereg$min_genes %||% 5
            tf_regulons <- do_call("GetTFRegulons", list(seurat_obj = srtobj, wgcna_name = wgcna_name))
            combined_output <- data.frame()
            for (cur_tf in unique(tf_regulons$tf)) {
                cur_targets <- do_call("GetTFTargetGenes", list(
                    seurat_obj = srtobj, selected_tfs = cur_tf, depth = depth,
                    use_regulons = use_regulons, wgcna_name = wgcna_name
                ))
                for (cur_type in c("positive", "negative")) {
                    cur_genes <- unique(subset(cur_targets, (Cor > 0) == (cur_type == "positive"))$gene)
                    if (length(cur_genes) < min_genes) next
                    enriched <- RunEnrichment(cur_genes, dbs = dbs, style = "enrichr")
                    if (is.null(enriched) || nrow(enriched) == 0) next
                    enriched$tf <- cur_tf
                    enriched$target_type <- cur_type
                    combined_output <- rbind(combined_output, enriched)
                }
            }
            if (nrow(combined_output) > 0) {
                srtobj <- do_call("SetEnrichrRegulonTable", list(
                    seurat_obj = srtobj, enrich_table = combined_output, wgcna_name = wgcna_name
                ))
            }
            enrich_df <- GetEnrichrRegulonTable(srtobj, wgcna_name = wgcna_name)
            if (!is.null(enrich_df) && nrow(enrich_df) > 0) {
                write_table(enrich_df, file.path(dreg_dir, "regulon_enrichment.tsv"))
                reporter$add2(
                    list(name = "Enrichment", contents = list(
                        list(kind = "descr", content = "The enrichment results of the TF regulons."),
                        list(kind = "table", src = file.path(dreg_dir, "regulon_enrichment.tsv"), data = list(nrows = 100))
                    )),
                    hs = c("TF Network", "Regulon Enrichment"),
                    ui = "tabs"
                )
            }
        }

        ############## ModuleRegulatoryHeatmap / ModuleRegulatoryNetworkPlot ##############
        # ModuleRegulatoryHeatmap uses fct_rev() in its aes, which resolves at
        # plot-build time in the hdWGCNA namespace; hdWGCNA does not import
        # forcats, so attach it here (the tutorial runs under tidyverse)
        suppressPackageStartupMessages(library(forcats))
        if (!is.null(tf_network$regulatory_heatmap)) {
            rhcase <- tf_network$regulatory_heatmap
            if (is.null(rhcase$feature)) stop("`feature` is required for ModuleRegulatoryHeatmap!")
            rhinfo <- case_info("Module Regulatory Heatmap::Regulatory Heatmap", tf_plots_dir, is_dir = FALSE, create = TRUE)
            devpars <- rhcase$devpars %||% list()
            more_formats <- rhcase$more_formats %||% list()
            descr <- rhcase$descr %||% "The module regulatory heatmap showing the TF regulatory activities on the modules."
            rhcase <- rhcase[setdiff(names(rhcase), c("devpars", "more_formats", "descr"))]
            p <- do_call("ModuleRegulatoryHeatmap", c(
                list(seurat_obj = srtobj),
                hd_args(c(list(wgcna_name = wgcna_name), rhcase), "ModuleRegulatoryHeatmap")
            ))
            save_plot(p, rhinfo$prefix, devpars, formats = c("png", more_formats))
            reporter$add2(
                list(name = "Plot", contents = list(
                    list(kind = "descr", content = descr),
                    reporter$image(rhinfo$prefix, more_formats, FALSE, kind = "image")
                )),
                hs = c("Module Regulatory Heatmap", "Regulatory Heatmap"),
                ui = "tabs"
            )
        }
        if (!is.null(tf_network$regulatory_network)) {
            rncase <- tf_network$regulatory_network
            if (is.null(rncase$feature)) stop("`feature` is required for ModuleRegulatoryNetworkPlot!")
            rninfo <- case_info("Module Regulatory Network::Regulatory Network", tf_plots_dir, is_dir = FALSE, create = TRUE)
            devpars <- rncase$devpars %||% list()
            more_formats <- rncase$more_formats %||% list()
            descr <- rncase$descr %||% "The module regulatory network plot showing the TF regulatory connections."
            rncase <- rncase[setdiff(names(rncase), c("devpars", "more_formats", "descr"))]
            p <- do_call("ModuleRegulatoryNetworkPlot", c(
                list(seurat_obj = srtobj),
                hd_args(c(list(wgcna_name = wgcna_name), rncase), "ModuleRegulatoryNetworkPlot")
            ))
            save_plot(p, rninfo$prefix, devpars, formats = c("png", more_formats))
            reporter$add2(
                list(name = "Plot", contents = list(
                    list(kind = "descr", content = descr),
                    reporter$image(rninfo$prefix, more_formats, FALSE, kind = "image")
                )),
                hs = c("Module Regulatory Network", "Regulatory Network"),
                ui = "tabs"
            )
        }
    }, error = function(e) {
        tb <- utils::tail(vapply(sys.calls(), function(x) paste(deparse(x), collapse = " "), character(1)), 15)
        log$warn(glue("Failed to run the TF network analysis: {conditionMessage(e)}"))
        # deparsed calls contain `{`/`}` which break the glue formatter of
        # poplog — strip them so a failing section never kills the job
        log$warn(glue("Traceback: {gsub('[{}]', '', paste(tb, collapse = ' || '))}"))
    })
}
