############## module preservation ##############
if (length(module_preservations) > 0) {
    if (is.null(ref_srtobj) && !ref_is_self) {
        log$warn("`module_preservations` requires `ref_srtobj` to be set! Skipped.")
    } else {
        # for a self-preservation test the reference is the query object
        # itself (with the network built by the core analysis)
        ref_obj <- if (ref_is_self) srtobj else ref_srtobj
        module_preservations <- lapply(module_preservations, function(x) list_update(module_preservations_defaults, x))
        pres_dir <- file.path(tables_dir, "preservation")
        dir.create(pres_dir, recursive = TRUE, showWarnings = FALSE)
        for (name in names(module_preservations)) {
            case <- module_preservations[[name]]
            tryCatch({
                log$info(glue("Running module preservation `{name}` ..."))
                project_modules <- case$project_modules %||% TRUE
                preserve <- case$preserve %||% TRUE
                plot <- case$plot %||% TRUE
                lollipop <- case$lollipop %||% FALSE
                netrep <- case$netrep
                devpars <- case$devpars %||% list()
                more_formats <- case$more_formats %||% list()
                descr <- case$descr %||% "The module preservation test of the modules in the query object against the reference object."
                project_args <- case$project_args %||% list()
                preserve_args <- case$preserve_args %||% list()
                case <- case[setdiff(names(case), c("project_modules", "preserve", "project_args", "preserve_args", "plot", "lollipop", "netrep", "devpars", "more_formats", "descr"))]
                info <- case_info(paste0("Module Preservation::", name), pres_dir, is_dir = FALSE, create = TRUE)
                if (isTRUE(project_modules)) {
                    log$info(glue("Projecting the modules to the reference object for `{name}` ..."))
                    srtobj <<- do_call("ProjectModules", c(
                        list(seurat_obj = srtobj, seurat_ref = ref_obj, wgcna_name = wgcna_name),
                        list_update(list(wgcna_name_proj = name), project_args)
                    ))
                }
                pres <- NULL
                if (isTRUE(preserve)) {
                    log$info(glue("Testing the module preservation for `{name}` ..."))
                    srtobj <<- do_call("ModulePreservation", c(
                        list(seurat_obj = srtobj, seurat_ref = ref_obj, wgcna_name = wgcna_name),
                        list_update(list(name = name), preserve_args)
                    ))
                    pres <- GetModulePreservation(srtobj, mod_name = name, wgcna_name = wgcna_name)
                    for (k in intersect(names(pres), c("Z", "obs", "p.values", "p.value"))) {
                        write_table(pres[[k]], paste0(info$prefix, ".", k, ".tsv"), row.names = TRUE)
                    }
                }
                report_items <- list(list(kind = "descr", content = descr))
                if (isTRUE(plot) && !is.null(pres)) {
                    ps <- do_call("PlotModulePreservation", c(
                        list(seurat_obj = srtobj, name = name, wgcna_name = wgcna_name),
                        hd_args(case, "PlotModulePreservation")
                    ))
                    save_plot_list(ps, info$prefix, auto_devpars(devpars, length(ps), ncol = 2), ncol = 2, formats = c("png", more_formats))
                    report_items <- c(report_items, list(reporter$image(info$prefix, more_formats, FALSE, kind = "image")))
                }
                if (isTRUE(lollipop) && !is.null(pres)) {
                    p <- do_call("PlotModulePreservationLollipop", c(
                        list(seurat_obj = srtobj, name = name, wgcna_name = wgcna_name),
                        hd_args(case, "PlotModulePreservationLollipop")
                    ))
                    save_plot(p, paste0(info$prefix, ".lollipop"), devpars, formats = c("png", more_formats))
                    report_items <- c(report_items, list(reporter$image(paste0(info$prefix, ".lollipop"), more_formats, FALSE, kind = "image")))
                }
                if (length(netrep) > 0) {
                    log$info(glue("Running the network reproduction test for `{name}` ..."))
                    netrep_args <- list_update(
                        list(name = name, wgcna_name = wgcna_name, wgcna_name_ref = wgcna_name),
                        netrep$args %||% list()
                    )
                    srtobj <<- do_call("ModulePreservationNetRep", c(
                        list(seurat_query = srtobj, seurat_ref = ref_obj),
                        netrep_args
                    ))
                    topo_dir <- file.path(plots_dir, "Topology")
                    dir.create(topo_dir, recursive = TRUE, showWarnings = FALSE)
                    for (tname in names(netrep$topology_heatmap %||% list())) {
                        tcase <- netrep$topology_heatmap[[tname]]
                        tryCatch({
                            if (is.null(tcase$mod)) stop("`mod` is required for the topology heatmap!")
                            use_ref <- isTRUE(tcase$use_ref)
                            tdevpars <- tcase$devpars %||% devpars
                            tmore <- tcase$more_formats %||% more_formats
                            tdescr <- tcase$descr %||% glue("The topology heatmap of module `{tcase$mod}`.")
                            tinfo <- case_info(paste0("Module Preservation::", paste0(name, " (", tname, ")")), topo_dir, is_dir = FALSE, create = TRUE)
                            tcase <- tcase[setdiff(names(tcase), c("use_ref", "descr", "devpars", "more_formats"))]
                            obj <- if (use_ref) ref_obj else srtobj
                            wn <- if (use_ref) netrep_args$wgcna_name_ref else wgcna_name
                            topo_args <- hd_args(tcase, "ModuleTopologyHeatmap")
                            if ("wgcna_name" %in% names(formals(get("ModuleTopologyHeatmap", asNamespace("hdWGCNA"))))) {
                                topo_args$wgcna_name <- wn
                            }
                            p <- do_call("ModuleTopologyHeatmap", c(list(seurat_obj = obj), topo_args))
                            save_plot(p, tinfo$prefix, tdevpars, formats = c("png", tmore))
                            reporter$add2(
                                list(name = "Plot", contents = list(
                                    list(kind = "descr", content = tdescr),
                                    reporter$image(tinfo$prefix, tmore, FALSE, kind = "image")
                                )),
                                hs = hs_section("Module Preservation", paste0(name, " (", tname, ")")),
                                ui = "tabs"
                            )
                        }, error = function(e2) {
                            log$warn(glue("Failed to plot topology heatmap `{tname}`: {conditionMessage(e2)}"))
                        })
                    }
                    for (tname in names(netrep$topology_barplot %||% list())) {
                        tcase <- netrep$topology_barplot[[tname]]
                        tryCatch({
                            if (is.null(tcase$mod)) stop("`mod` is required for the topology barplot!")
                            use_ref <- isTRUE(tcase$use_ref)
                            tdevpars <- tcase$devpars %||% devpars
                            tmore <- tcase$more_formats %||% more_formats
                            tdescr <- tcase$descr %||% glue("The topology barplot of module `{tcase$mod}`.")
                            tinfo <- case_info(paste0("Module Preservation::", paste0(name, " (", tname, ")")), topo_dir, is_dir = FALSE, create = TRUE)
                            tcase <- tcase[setdiff(names(tcase), c("use_ref", "descr", "devpars", "more_formats"))]
                            obj <- if (use_ref) ref_obj else srtobj
                            wn <- if (use_ref) netrep_args$wgcna_name_ref else wgcna_name
                            topo_args <- hd_args(tcase, "ModuleTopologyBarplot")
                            if ("wgcna_name" %in% names(formals(get("ModuleTopologyBarplot", asNamespace("hdWGCNA"))))) {
                                topo_args$wgcna_name <- wn
                            }
                            p <- do_call("ModuleTopologyBarplot", c(list(seurat_obj = obj), topo_args))
                            save_plot(p, tinfo$prefix, tdevpars, formats = c("png", tmore))
                            reporter$add2(
                                list(name = "Plot", contents = list(
                                    list(kind = "descr", content = tdescr),
                                    reporter$image(tinfo$prefix, tmore, FALSE, kind = "image")
                                )),
                                hs = hs_section("Module Preservation", paste0(name, " (", tname, ")")),
                                ui = "tabs"
                            )
                        }, error = function(e2) {
                            log$warn(glue("Failed to plot topology barplot `{tname}`: {conditionMessage(e2)}"))
                        })
                    }
                }
                table_items <- list()
                if (!is.null(pres)) {
                    table_items <- list(
                        list(kind = "descr", content = "The module preservation statistics."),
                        list(kind = "table", src = paste0(info$prefix, ".Z.tsv"), data = list(nrows = 100)),
                        list(kind = "table", src = paste0(info$prefix, ".obs.tsv"), data = list(nrows = 100))
                    )
                }
                # only add the case-level section when it has actual content;
                # a pure-netrep case (no preserve/plot) would otherwise add
                # an empty tab (its topology plots are reported separately)
                tabs <- list()
                if (length(report_items) > 1) {
                    tabs <- c(tabs, list(list(name = "Plot", contents = report_items)))
                }
                if (length(table_items) > 0) {
                    tabs <- c(tabs, list(list(name = "Tables", contents = table_items)))
                }
                if (length(tabs) > 0) {
                    # `add2` takes one tab per argument (`...`); a list of
                    # tabs passed as a single argument would be nested one
                    # level too deep and break the report render
                    do.call(reporter$add2, c(tabs, list(hs = hs_section("Module Preservation", name), ui = "tabs")))
                }
            }, error = function(e) {
                log$warn(glue("Failed to run module preservation `{name}`: {conditionMessage(e)}"))
            })
        }
    }
}
