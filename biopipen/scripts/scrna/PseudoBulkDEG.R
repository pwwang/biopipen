library(rlang)
library(dplyr)
library(plotthis)
library(biopipen.utils)

sobjfile <- {{in.sobjfile | r}}
outdir <- {{out.outdir | r}}
joboutdir <- {{job.outdir | r}}
each <- {{envs.each | r}}
subset <- {{envs.subset | r}}
ncores <- {{envs.ncores | r}}
mutaters <- {{envs.mutaters | r}}
cache <- {{ envs.cache | r }}
aggregate_by <- {{envs.aggregate_by | r}}
layer <- {{envs.layer | r}}
assay <- {{envs.assay | r}}
error <- {{ envs.error | r }}
group_by <- {{envs.group_by | default: envs[["group-by"]] | default: None | r}}
ident_1 <- {{envs.ident_1 | default: envs[["ident-1"]] | default: None | r}}
ident_2 <- {{envs.ident_2 | default: envs[["ident-2"]] | default: None | r}}
each <- {{ envs.each | r }}
dbs <- {{ envs.dbs | r }}
sigmarkers <- {{ envs.sigmarkers | r }}
enrich_style <- {{ envs.enrich_style | r }}
paired_by <- {{envs.paired_by | default: envs[["paired-by"]] | default: None | r}}
tool <- {{envs.tool | r}}
allmarker_plots_defaults <- {{ envs.allmarker_plots_defaults | r }}
allmarker_plots <- {{ envs.allmarker_plots | r }}
allenrich_plots_defaults <- {{ envs.allenrich_plots_defaults | r }}
allenrich_plots <- {{ envs.allenrich_plots | r }}
marker_plots_defaults <- {{ envs.marker_plots_defaults | r }}
marker_plots <- {{ envs.marker_plots | r }}
enrich_plots_defaults <- {{ envs.enrich_plots_defaults | r }}
enrich_plots <- {{ envs.enrich_plots | r }}
overlaps_defaults <- {{ envs.overlaps_defaults | r }}
overlaps <- {{ envs.overlaps | r }}
cases <- {{envs.cases | r}}

aggregate_by <- unique(c(aggregate_by, group_by, paired_by, each))
if (isTRUE(cache)) { cache <- joboutdir }

log <- get_logger()
reporter <- get_reporter()

log$info("Loading Seurat object ...")
srtobj <- read_obj(sobjfile)

if (!is.null(mutaters) && length(mutaters) > 0) {
    log$info("Mutating metadata columns ...")
    srtobj@meta.data <- srtobj@meta.data %>% mutate(!!!lapply(mutaters, parse_expr))
}

allmarker_plots <- lapply(allmarker_plots, function(x) {
    list_update(allmarker_plots_defaults, x)
})
allenrich_plots <- lapply(allenrich_plots, function(x) {
    list_update(allenrich_plots_defaults, x)
})
marker_plots <- lapply(marker_plots, function(x) {
    list_update(marker_plots_defaults, x)
})
enrich_plots <- lapply(enrich_plots, function(x) {
    list_update(enrich_plots_defaults, x)
})
overlaps <- lapply(overlaps, function(x) {
    list_update(overlaps_defaults, x)
})

defaults <- list(
    each = each,
    error = error,
    subset = subset,
    aggregate_by = aggregate_by,
    layer = layer,
    assay = assay %||% DefaultAssay(srtobj),
    group_by = group_by,
    ident_1 = ident_1,
    ident_2 = ident_2,
    dbs = dbs,
    ncores = ncores,
    sigmarkers = sigmarkers,
    enrich_style = enrich_style,
    paired_by = paired_by,
    tool = tool,
    cache = cache,
    allmarker_plots_defaults = allmarker_plots_defaults,
    allmarker_plots = allmarker_plots,
    allenrich_plots_defaults = allenrich_plots_defaults,
    allenrich_plots = allenrich_plots,
    marker_plots_defaults = marker_plots_defaults,
    marker_plots = marker_plots,
    enrich_plots_defaults = enrich_plots_defaults,
    enrich_plots = enrich_plots,
    overlaps_defaults = overlaps_defaults,
    overlaps = overlaps
)

expand_each <- function(name, case) {
    outcases <- list()

    if (is.null(case$each) || is.na(case$each) || nchar(case$each) == 0 || isFALSE(each)) {
        if (length(case$allmarker_plots) > 0 || length(allenrich_plots) > 0 || length(overlaps) > 0) {
            stop("Cannot perform allmarker_plots, allenrich_plots, or overlaps without 'each' defined.")
        }

        case$aggregate_by <- unique(c(case$aggregate_by, case$group_by, case$paired_by))

        case$marker_plots <- lapply(
            case$marker_plots,
            function(x) { list_update(case$marker_plots_defaults, x) }
        )
        case$marker_plots_defaults <- NULL

        case$enrich_plots <- lapply(
            case$enrich_plots,
            function(x) { list_update(case$enrich_plots_defaults, x) }
        )
        case$enrich_plots_defaults <- NULL

        case$overlaps <- lapply(
            case$overlaps,
            function(x) { list_update(case$overlaps_defaults, x) }
        )
        case$overlaps_defaults <- NULL

        outcases[[name]] <- case
    } else {
        eachs <- if (!is.null(case$subset)) {
            srtobj@meta.data %>%
                filter(!!parse_expr(case$subset)) %>%
                pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
        } else {
            srtobj@meta.data %>%
                pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
        }

        if (length(cases) == 0 && name == "DEG Analysis") {
            name <- case$each
        } else {
            name <- paste0(name, " (", case$each, ")")
        }

        case$aggregate_by <- unique(c(case$aggregate_by, case$group_by, case$paired_by, case$each))

        for (each in eachs) {
            newname <- paste0(name, "::", each)
            newcase <- case

            newcase$original_case <- name
            newcase$each_name <- case$each
            newcase$each <- each

            newcase$alleach_plots_defaults <- NULL
            newcase$alleach_plots <- NULL
            newcase$original_subset <- case$subset

            if (!is.null(case$subset)) {
                newcase$subset <- paste0(case$subset, " & ", bQuote(case$each), " == '", each, "'")
            } else {
                newcase$subset <- paste0(bQuote(case$each), " == '", each, "'")
            }

            newcase$marker_plots <- lapply(
                case$marker_plots,
                function(x) { list_update(case$marker_plots_defaults, x) }
            )
            newcase$marker_plots_defaults <- NULL

            newcase$enrich_plots <- lapply(
                case$enrich_plots,
                function(x) { list_update(case$enrich_plots_defaults, x) }
            )
            newcase$enrich_plots_defaults <- NULL

            # Will be processed by the case itself, which collects the markers
            newcase$allmarker_plots <- NULL
            newcase$allmarker_plots_defaults <- NULL
            newcase$allenrich_plots <- NULL
            newcase$allenrich_plots_defaults <- NULL
            newcase$overlaps <- NULL
            newcase$overlaps_defaults <- NULL

            outcases[[newname]] <- newcase
        }

        if (length(case$overlaps) > 0 || length(case$allmarker_plots) > 0 || length(case$allenrich_plots) > 0) {
            ovcase <- case

            ovcase$allexprs <- list()
            ovcase$markers <- list()
            ovcase$allmarker_plots <- lapply(
                ovcase$allmarker_plots,
                function(x) { list_update(ovcase$allmarker_plots_defaults, x) }
            )
            ovcase$allmarker_plots_defaults <- NULL

            ovcase$enriches <- list()
            ovcase$allenrich_plots <- lapply(
                ovcase$allenrich_plots,
                function(x) { list_update(ovcase$allenrich_plots_defaults, x) }
            )
            ovcase$allenrich_plots_defaults <- NULL

            ovcase$overlaps <- lapply(
                ovcase$overlaps,
                function(x) { list_update(ovcase$overlaps_defaults, x) }
            )
            ovcase$overlaps_defaults <- NULL
            outcases[[name]] <- ovcase
        }
    }
    outcases
}
cases <- expand_cases(cases, defaults, expand_each, default_case = "DEG Analysis")

log$info("Running cases ...")

process_markers <- function(markers, info, case) {
    ## Attributes lost
    # markers <- markers %>%
    #     mutate(gene = as.character(gene)) %>%
    #     arrange(p_val_adj, desc(abs(avg_log2FC)))

    empty <- if (case$enrich_style == "enrichr") {
        data.frame(
            Database = character(0),
            Term = character(0),
            Overlap = character(0),
            P.value = numeric(0),
            Adjusted.P.value = numeric(0),
            Odds.Ratio = numeric(0),
            Combined.Score = numeric(0),
            Genes = character(0),
            Rank = numeric(0)
        )
    } else {  # clusterProfiler
        data.frame(
            ID = character(0),
            Description = character(0),
            GeneRatio = character(0),
            BgRatio = character(0),
            Count = integer(0),
            pvalue = numeric(0),
            p.adjust = numeric(0),
            qvalue = numeric(0),
            geneID = character(0),
            Database = character(0)
        )
    }
    if (is.null(markers) || nrow(markers) == 0) {
        if (case$error) {
            stop("Error: No markers found in case '", info$name, "'.")
        } else {
            log$warn("! Warning: No markers found in case '", info$name, "'.")
            reporter$add2(
                list(
                    name = "Warning",
                    contents = list(list(kind = "error", content = "No markers found.", kind_ = "warning"))),
                hs = c(info$section, info$name),
                hs2 = "DEG Analysis",
                ui = "tabs"
            )
            return(empty)
        }
    }
    markers$gene <- as.character(markers$gene)
    markers$p_val_adj <- as.numeric(markers$p_val_adj)
    markers$log2FC <- as.numeric(markers$log2FC)
    markers <- markers[order(markers$p_val_adj, -abs(markers$log2FC)), ]

    # Save markers
    write.table(markers, file.path(info$prefix, "degs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    sigmarkers <- markers %>% filter(!!parse_expr(case$sigmarkers))
    write.table(sigmarkers, file.path(info$prefix, "sigdegs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    reporter$add2(
        list(
            name = "Table",
            contents = list(
                list(kind = "descr", content = paste0(
                    "Showing top 100 DEGs ordered by p_val_adj ascendingly, then abs(log2FC) descendingly. ",
                    "Use 'Download the entire data' button to download all significant markers by '",
                    html_escape(case$sigmarkers), "'."
                )),
                list(kind = "table", src = file.path(info$prefix, "sigdegs.tsv"), data = list(nrows = 100))
            )
        ),
        hs = c(info$section, info$name),
        hs2 = ifelse(is.null(case$ident), "DEGs", paste0("DEGs (", case$ident, ")")),
        ui = "tabs"
    )

    for (plotname in names(case$marker_plots)) {
        plotargs <- case$marker_plots[[plotname]]
        plotargs$degs <- markers
        rownames(plotargs$degs) <- make.unique(markers$gene)
        plotargs$outprefix <- file.path(info$prefix, paste0("degs.", slugify(plotname)))
        do_call(VizBulkDEGs, plotargs)
        reporter$add2(
            list(
                name = plotname,
                contents = list(reporter$image(plotargs$outprefix, plotargs$more_formats, plotargs$save_code))),
            hs = c(info$section, info$name),
            hs2 = ifelse(is.null(case$ident), "DEGs", paste0("DEGs (", case$ident, ")")),
            ui = "tabs"
        )
    }

    # Do enrichment analysis
    significant_markers <- unique(sigmarkers$gene)
    empty <- if (case$enrich_style == "enrichr") {
        data.frame(
            Database = character(0),
            Term = character(0),
            Overlap = character(0),
            P.value = numeric(0),
            Adjusted.P.value = numeric(0),
            Odds.Ratio = numeric(0),
            Combined.Score = numeric(0),
            Genes = character(0),
            Rank = numeric(0)
        )
    } else {  # clusterProfiler
        data.frame(
            ID = character(0),
            Description = character(0),
            GeneRatio = character(0),
            BgRatio = character(0),
            Count = integer(0),
            pvalue = numeric(0),
            p.adjust = numeric(0),
            qvalue = numeric(0),
            geneID = character(0),
            Database = character(0)
        )
    }

    if (length(significant_markers) < 5) {
        if (case$error) {
            stop("Error: Not enough significant DEGs with '", case$sigmarkers, "' in case '", info$name, "' found (< 5) for enrichment analysis.")
        } else {
            message <- paste0("Not enough significant DEGs with '", case$sigmarkers, "' found (< 5) for enrichment analysis.")
            log$warn("! Error: {message}")
            reporter$add2(
                list(
                    name = "Warning",
                    contents = list(list(kind = "error", content = message, kind_ = "warning"))),
                hs = c(info$section, info$name),
                hs2 = "Enrichment Analysis",
                ui = "tabs"
            )
        }
        return(empty)
    } else {
        tryCatch({
            enrich <- RunEnrichment(
                significant_markers,
                dbs = case$dbs, style = case$enrich_style)

            write.table(enrich, file.path(info$prefix, "enrich.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
            reporter$add2(
                list(
                    name = "Table",
                    contents = list(list(kind = "table", src = file.path(info$prefix, "enrich.tsv"), data = list(nrows = 100)))
                ),
                hs = c(info$section, info$name),
                hs2 = "Enrichment Analysis",
                ui = "tabs"
            )

            # Visualize enriched terms
            if (length(case$enrich_plots) > 0) {
                for (db in case$dbs) {
                    plots <- list()
                    for (plotname in names(case$enrich_plots)) {
                        plotargs <- case$enrich_plots[[plotname]]
                        plotargs$data <- enrich[enrich$Database == db, , drop = FALSE]

                        p <- do_call(VizEnrichment, plotargs)

                        if (plotargs$plot_type == "bar") {
                            attr(p, "height") <- attr(p, "height") / 1.5
                        }
                        outprefix <- file.path(info$prefix, paste0("enrich.", slugify(db), ".", slugify(plotname)))
                        save_plot(p, outprefix, plotargs$devpars, formats = "png")
                        plots[[length(plots) + 1]] <- reporter$image(outprefix, c(), FALSE)
                    }
                    reporter$add2(
                        list(name = db, contents = plots),
                        hs = c(info$section, info$name),
                        hs2 = "Enrichment Analysis",
                        ui = "tabs"
                    )
                }
            }
            return(enrich)
        }, error = function(e) {
            if (case$error) {
                stop("Error: ", e$message)
            } else {
                log$warn("! Error: {e$message}")
                reporter$add2(
                    list(
                        name = "Warning",
                        contents = list(list(kind = "error", content = e$message, kind_ = "warning"))),
                    hs = c(info$section, info$name),
                    hs2 = "Enrichment Analysis",
                    ui = "tabs"
                )
            }
            return(empty)
        })
    }
}

process_allmarkers <- function(markers, plotcases, casename, groupname) {
    name <- paste0(casename, "::", paste0(groupname, " (All DEGs)"))
    info <- case_info(name, outdir, create = TRUE)

    for (plotname in names(plotcases)) {
        plotargs <- plotcases[[plotname]]
        plotargs$degs <- markers
        plotargs$outprefix <- file.path(info$prefix, slugify(plotname))
        do_call(VizBulkDEGs, plotargs)
        reporter$add2(
            list(
                name = plotname,
                contents = list(reporter$image(plotargs$outprefix, plotargs$more_formats, plotargs$save_code))
            ),
            hs = c(info$section, info$name),
            ui = "tabs"
        )
    }
}

process_allenriches <- function(enriches, plotcases, casename, groupname) {
    name <- paste0(casename, "::", paste0(groupname, " (All Enrichments)"))
    info <- case_info(name, outdir, create = TRUE)
    dbs <- unique(as.character(enriches$Database))

    for (db in dbs) {
        plots <- list()
        for (plotname in names(plotcases)) {
            plotargs <- plotcases[[plotname]]
            plotargs <- extract_vars(plotargs, "devpars")
            plotargs$data <- enriches[enriches$Database == db, , drop = FALSE]
            if (plotargs$plot_type == "heatmap") {
                plotargs$group_by <- groupname
                plotargs$show_row_names = plotargs$show_row_names %||% TRUE
                plotargs$show_column_names = plotargs$show_column_names %||% TRUE
            }

            p <- do_call(VizEnrichment, plotargs)

            if (plotargs$plot_type == "bar") {
                attr(p, "height") <- attr(p, "height") / 1.5
            }
            outprefix <- file.path(info$prefix, paste0("allenrich.", slugify(db), ".", slugify(plotname)))
            save_plot(p, outprefix, devpars, formats = "png")
            plots[[length(plots) + 1]] <- reporter$image(outprefix, c(), FALSE)
        }
        reporter$add2(
            list(name = db, contents = plots),
            hs = c(info$section, info$name),
            hs2 = plotname,
            ui = "tabs"
        )
    }
}

process_overlaps <- function(markers, ovcases, casename, groupname) {
    name <- paste0(casename, "::", paste0(groupname, " (Overlaps)"))
    info <- case_info(name, outdir, create = TRUE)

    for (plotname in names(ovcases)) {
        args <- extract_vars(
            ovcases[[plotname]],
            sigm = "sigmarkers", "more_formats", "save_code", "devpars", "plot_type",
            allow_nonexisting = TRUE
        )

        sigm <- sigm %||% sigmarkers
        ugroups <- unique(markers[[groupname]])
        m <- lapply(ugroups, function(g) {
            markers[markers[[groupname]] == g, , drop = FALSE] %>%
                filter(!!parse_expr(sigm)) %>%
                pull("gene") %>% unique()
        })
        names(m) <- ugroups

        if (plot_type == "venn") {
            args$data <- m
            args$in_form <- "list"
            prefix <- file.path(info$prefix, slugify(plotname))
            p <- do_call(gglogger::register(VennDiagram), args)
            save_plot(p, prefix, devpars, formats = c("png", more_formats))
            if (save_code) {
                save_plotcode(
                    p, prefix,
                    c("library(plotthis)", "load('data.RData')", "invisible(list2env(args, .GlobalEnv))"),
                    "args",
                    auto_data_setup = FALSE
                )
            }
        } else {
            args$data <- m
            args$in_form <- "list"
            prefix <- file.path(info$prefix, slugify(plotname))
            p <- do_call(gglogger::register(UpsetPlot), args)
            save_plot(p, prefix, devpars, formats = c("png", more_formats))
            if (save_code) {
                save_plotcode(
                    p, prefix,
                    c("library(plotthis)", "load('data.RData')", "invisible(list2env(args, .GlobalEnv))"),
                    "args",
                    auto_data_setup = FALSE
                )
            }
        }

        reporter$add2(
            list(
                name = plotname,
                contents = list(reporter$image(prefix, more_formats, save_code))
            ),
            hs = c(info$section, info$name),
            ui = "tabs"
        )
    }
}

run_case <- function(name) {
    case <- cases[[name]]
    log$info("----------------------------------------")
    log$info("Case: {name} ...")

    case <- extract_vars(
        case,
        "dbs", "sigmarkers", "allmarker_plots", "allenrich_plots", "marker_plots", "enrich_plots",
        "overlaps", "original_case", "markers", "enriches", "each_name", "each", "enrich_style",
        "aggregate_by", "subset", "layer", "assay", "group_by", "ident_1", "ident_2", "original_subset",
        "paired_by", "tool", "error", "ncores", "cache", "allexprs",
        allow_nonexisting = TRUE
    )

    if (!is.null(markers) || !is.null(enriches)) {
        if (!is.null(markers) && length(markers) > 0) {
            log$info("Summarizing DEGs in subcases (by each: {each}) ...")
            # handle the overlaps / allmarkers analysis here
            if (!is.data.frame(markers)) {
                each_levels <- names(markers)
                markers <- do_call(rbind, lapply(each_levels, function(x) {
                    markers_df <- markers[[x]]
                    if (is.null(markers_df) || nrow(markers_df) == 0) {
                        return(NULL)
                    }
                    if (nrow(markers_df) > 0) {
                        markers_df[[each]] <- x
                    } else {
                        markers_df[[each]] <- character(0)  # Empty case
                    }
                    markers_df
                }))
                markers[[each]] <- factor(markers[[each]], levels = each_levels)
            }
            # gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, diff_pct, <each>

            if (!is.data.frame(allexprs)) {
                meta <- do_call(rbind, lapply(allexprs, attr, "meta"))
                allexprs <- do_call(cbind, allexprs)
            } else {
                meta <- attr(allexprs, "meta")
            }

            if (length(allmarker_plots) > 0) {
                log$info("Visualizing all DEGs together ...")
                attr(markers, "object") <- allexprs
                attr(markers, "meta") <- meta
                attr(markers, "group_by") <- each
                attr(markers, "paired_by") <- paired_by
                attr(markers, "ident_1") <- NULL
                attr(markers, "ident_2") <- NULL
                process_allmarkers(markers, allmarker_plots, name, each)
            }

            if (length(overlaps) > 0) {
                log$info("Visualizing overlaps between subcases ...")
                process_overlaps(markers, overlaps, name, each)
            }

        }

        if (!is.null(enriches) && length(enriches) > 0) {
            log$info("Summarizing enrichments in subcases (by each: {each}) ...")
            if (!is.data.frame(enriches)) {
                each_levels <- names(enriches)
                enriches <- do_call(rbind, lapply(each_levels, function(x) {
                    enrich_df <- enriches[[x]]
                    if (is.null(enrich_df) || nrow(enrich_df) == 0) {
                        return(NULL)
                    }
                    if (nrow(enrich_df) > 0) {
                        enrich_df[[each]] <- x
                    } else {
                        enrich_df[[each]] <- character(0)  # Empty case
                    }
                    enrich_df
                }))
                if (!is.null(enriches) && nrow(enriches) > 0) {
                    enriches[[each]] <- factor(enriches[[each]], levels = each_levels)
                }
            }

            if (length(allenrich_plots) > 0 && !is.null(enriches) && nrow(enriches) > 0) {
                log$info("Visualizing all enrichments together ...")
                process_allenriches(enriches, allenrich_plots, name, each)
            }
        }

        return(invisible())
    }

    info <- case_info(name, outdir, create = TRUE)
    exprs <- AggregateExpressionPseudobulk(
        srtobj, aggregate_by = aggregate_by, layer = layer, assay = assay,
        subset = subset, log = log
    )
    markers <- tryCatch(
        {
            RunDEGAnalysis(
                exprs, group_by = group_by, ident_1 = ident_1, ident_2 = ident_2,
                paired_by = paired_by, tool = tool, log = log, ncores = ncores,
                cache = cache
            )
        }, error = function(e) {
            if (error) {
                stop("Error: ", e$message)
            } else {
                log$warn("! Error: {e$message}")
                reporter$add2(
                    list(
                        name = "Warning",
                        contents = list(list(kind = "error", content = e$message, kind_ = "warning"))),
                    hs = c(info$section, info$name),
                    hs2 = "DEG Analysis",
                    ui = "tabs"
                )
                return(invisible())
            }
        }
    )

    enrich <- process_markers(markers, info = info, case = list(
        dbs = dbs,
        sigmarkers = sigmarkers,
        enrich_style = enrich_style,
        marker_plots = marker_plots,
        enrich_plots = enrich_plots,
        error = error,
        ident = if (is.null(case$ident_2)) case$ident_1 else paste0(case$ident_1, " vs ", case$ident_2)
    ))

    if (!is.null(original_case) && !is.null(cases[[original_case]])) {
        if (!is.null(markers)) {
            markers[[each_name]] <- each
        }
        cases[[original_case]]$markers[[each]] <<- markers
        cases[[original_case]]$enriches[[each]] <<- enrich
        cases[[original_case]]$allexprs[[each]] <<- exprs
    }

    invisible()
}

sapply(names(cases), run_case)

reporter$save(joboutdir)
