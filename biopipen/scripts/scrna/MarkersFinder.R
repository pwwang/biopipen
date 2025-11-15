library(rlang)
library(dplyr)
library(Seurat)
library(tidyseurat)
library(plotthis)
library(biopipen.utils)

log <- get_logger()
reporter <- get_reporter()

srtfile <- {{ in.srtobj | r }}
outdir <- {{ out.outdir | r }}
joboutdir <- {{ job.outdir | r }}

ncores <- {{ envs.ncores | int }}
mutaters <- {{ envs.mutaters | r }}
group_by <- {{ envs.group_by | default: envs["group-by"] | default: None | r }}
ident_1 <- {{ envs.ident_1 | default: envs["ident-1"] | default: None | r }}
ident_2 <- {{ envs.ident_2 | default: envs["ident-2"] | default: None | r }}
each <- {{ envs.each | r }}
dbs <- {{ envs.dbs | r }}
sigmarkers <- {{ envs.sigmarkers | r }}
enrich_style <- {{ envs.enrich_style | r }}
assay <- {{ envs.assay | r }}
error <- {{ envs.error | r }}
subset <- {{ envs.subset | r }}
cache <- {{ envs.cache | r }}
rest <- {{ envs.rest | r: todot="-" }}
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
cases <- {{ envs.cases | r: todot="-", skip=1 }}

if (isTRUE(cache)) { cache <- joboutdir }

set.seed(8525)
if (ncores > 1) {
    options(future.globals.maxSize = Inf)
    plan(strategy = "multicore", workers = ncores)
}

log$info("Reading Seurat object ...")
srtobj <- read_obj(srtfile)


if (!is.null(mutaters) && length(mutaters) > 0) {
    log$info("Mutating meta data ...")
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
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
    group_by = group_by,
    ident_1 = ident_1,
    ident_2 = ident_2,
    dbs = dbs,
    sigmarkers = sigmarkers,
    enrich_style = enrich_style,
    assay = assay %||% DefaultAssay(srtobj),
    each = each,
    error = error,
    subset = subset,
    allmarker_plots_defaults = allmarker_plots_defaults,
    allmarker_plots = allmarker_plots,
    allenrich_plots_defaults = allenrich_plots_defaults,
    allenrich_plots = allenrich_plots,
    marker_plots_defaults = marker_plots_defaults,
    marker_plots = marker_plots,
    enrich_plots_defaults = enrich_plots_defaults,
    enrich_plots = enrich_plots,
    overlaps_defaults = overlaps_defaults,
    overlaps = overlaps,
    cache = cache,
    rest = rest
)

log$info("Expanding cases ...")

post_casing <- function(name, case) {
    outcases <- list()

    case$group_by <- case$group_by %||% GetIdentityColumn(srtobj)

    if (is.null(case$each) || is.na(case$each) || nchar(case$each) == 0 || isFALSE(each)) {
        # single cases, no need to expand
        if (length(case$ident_1) > 0 && length(case$overlaps) > 0) {
            stop("Cannot perform 'overlaps' with a single comparison (ident-1 is set) in case '", name, "'")
        }
        if (length(case$ident_1) > 0 && length(case$allmarker_plots) > 0) {
            stop("Cannot perform 'allmarker_plots' with a single comparison (ident-1 is set) in case '", name, "'")
        }
        if (length(case$ident_1) > 0 && length(case$allenrich_plots) > 0) {
            stop("Cannot perform 'allenrich_plots' with a single comparison (ident-1 is set) in case '", name, "'")
        }

        case$allmarker_plots <- lapply(
            case$allmarker_plots,
            function(x) { list_update(case$allmarker_plots_defaults, x) }
        )
        case$allmarker_plots_defaults <- NULL

        case$allenrich_plots <- lapply(
            case$allenrich_plots,
            function(x) { list_update(case$allenrich_plots_defaults, x) }
        )
        case$allenrich_plots_defaults <- NULL

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
    } else {  # !no_each
        eachs <- if (!is.null(case$subset)) {
            srtobj@meta.data %>%
                filter(!!parse_expr(case$subset)) %>%
                pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
        } else {
            srtobj@meta.data %>%
                pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
        }
        if (length(case$overlaps) > 0 && is.null(case$ident_1)) {
            stop("Cannot perform 'overlaps' analysis with 'each' and without 'ident_1' in case '", name, "'")
        }

        if (length(cases) == 0 && name == "Marker Discovery") {
            name <- case$each
        } else {
            name <- paste0(name, " (", case$each, ")")
        }

        for (each in eachs) {
            newname <- paste0(name, "::", each)
            newcase <- case

            newcase$original_case <- name
            newcase$each_name <- case$each
            newcase$each <- each
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
cases <- expand_cases(cases, defaults, post_casing, default_case = "Marker Discovery")

log$info("Running cases ...")

process_markers <- function(markers, info, case) {
    ## Attributes lost
    # markers <- markers %>%
    #     mutate(gene = as.character(gene)) %>%
    #     arrange(p_val_adj, desc(abs(avg_log2FC)))
    markers$gene <- as.character(markers$gene)
    markers <- markers[order(markers$p_val_adj, -abs(markers$avg_log2FC)), ]

    # Save markers
    write.table(markers, file.path(info$prefix, "markers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    sigmarkers <- markers %>% filter(!!parse_expr(case$sigmarkers))
    write.table(sigmarkers, file.path(info$prefix, "sigmarkers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    reporter$add2(
        list(
            name = "Table",
            contents = list(
                list(kind = "descr", content = paste0(
                    "Showing top 100 markers ordered by p_val_adj ascendingly, then abs(avg_log2FC) descendingly. ",
                    "Use 'Download the entire data' button to download all significant markers by '",
                    html_escape(case$sigmarkers), "'."
                )),
                list(kind = "table", src = file.path(info$prefix, "sigmarkers.tsv"), data = list(nrows = 100))
            )
        ),
        hs = c(info$section, info$name),
        hs2 = ifelse(is.null(case$ident), "Markers", paste0("Markers (", case$ident, ")")),
        ui = "tabs"
    )

    if (nrow(markers) > 0) {
        for (plotname in names(case$marker_plots)) {
            plotargs <- case$marker_plots[[plotname]]
            plotargs$markers <- markers
            plotargs$object <- case$object
            plotargs$comparison_by <- case$group_by
            plotargs$outprefix <- file.path(info$prefix, paste0("markers.", slugify(plotname)))
            do_call(VizDEGs, plotargs)
            reporter$add2(
                list(
                    name = plotname,
                    contents = list(reporter$image(plotargs$outprefix, plotargs$more_formats, plotargs$save_code))),
                hs = c(info$section, info$name),
                hs2 = ifelse(is.null(case$ident), "Markers", paste0("Markers (", case$ident, ")")),
                ui = "tabs"
            )
        }
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
            stop("Error: Not enough significant markers with '", case$sigmarkers, "' in case '", info$name, "' found (< 5) for enrichment analysis.")
        } else {
            message <- paste0("Not enough significant markers with '", case$sigmarkers, "' found (< 5) for enrichment analysis.")
            log$warn("  ! Error: {message}")
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
                        plotargs <- extract_vars(case$enrich_plots[[plotname]], "descr", allow_nonexisting = TRUE)
                        plotargs$data <- enrich[enrich$Database == db, , drop = FALSE]

                        p <- tryCatch(
                            do_call(VizEnrichment, plotargs),
                            error = function(e) {
                                stop("Failed to plot enrichment for database '", db, "' with plot '", plotname, "': ", e$message)
                            }
                        )

                        if (plotargs$plot_type == "bar") {
                            attr(p, "height") <- attr(p, "height") / 1.5
                            descr <- descr %||% glue::glue(
                                "The bar plot shows the top enriched terms in database '{db}', ",
                                "the x-axis shows the -log10 of the adjusted p-values, ",
                                "and the y-axis shows the term names. The number next to each bar indicates the overlap gene count."
                            )
                        }
                        outprefix <- file.path(info$prefix, paste0("enrich.", slugify(db), ".", slugify(plotname)))
                        save_plot(p, outprefix, plotargs$devpars, formats = "png")
                        if (!is.null(descr)) {
                            plots[[length(plots) + 1]] <- list(kind = "descr", content = glue::glue(descr))
                        }
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
                log$warn("  ! Error: {e$message}")
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

process_allmarkers <- function(markers, object, comparison_by, plotcases, casename, groupname, subset_by_group = TRUE) {
    name <- paste0(casename, "::", paste0(groupname, " (All Markers)"))
    info <- case_info(name, outdir, create = TRUE)

    for (plotname in names(plotcases)) {
        log$info("  {plotname} ...")
        plotargs <- plotcases[[plotname]]
        plotargs$markers <- markers
        plotargs$object <- object
        plotargs$comparison_by <- comparison_by
        if (subset_by_group)
            plotargs$subset_by <- groupname
        plotargs$outprefix <- file.path(info$prefix, slugify(plotname))
        do_call(VizDEGs, plotargs)
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
            log$info("  {plotname} ({db}) ...")
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
        log$info("  {plotname} ...")
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
    log$info("Case: {name} ...")

    case <- extract_vars(
        case,
        "dbs", "sigmarkers", "allmarker_plots", "allenrich_plots", "marker_plots", "enrich_plots",
        "overlaps", "original_case", "markers", "enriches", "each_name", "each", "enrich_style", "original_subset",
        subset_ = "subset",
        allow_nonexisting = TRUE
    )

    if (!is.null(markers) || !is.null(enriches)) {
        if (!is.null(markers)) {  # It is the overlap/allmarker case
            log$info("- Summarizing markers in subcases (by each: {each}) ...")
            # handle the overlaps / allmarkers analysis here
            if (!is.data.frame(markers)) {
                each_levels <- names(markers)
                markers <- do_call(rbind, lapply(each_levels, function(x) {
                    markers_df <- markers[[x]]
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

            if (length(allmarker_plots) > 0) {
                log$info("- Visualizing all markers together ...")
                if (is.null(original_subset)) {
                    attr(markers, "object") <- srtobj
                } else {
                    attr(markers, "object") <- filter(srtobj, !!parse_expr(original_subset))
                }
                attr(markers, "group_by") <- each
                attr(markers, "ident_1") <- NULL
                attr(markers, "ident_2") <- NULL
                if (!is.null(markers) && nrow(markers) > 0) {
                    process_allmarkers(
                        markers,
                        object = if (is.null(original_subset)) srtobj else filter(srtobj, !!parse_expr(original_subset)),
                        comparison_by = group_by,
                        allmarker_plots,
                        name,
                        each
                    )
                }
            }

            if (length(overlaps) > 0) {
                log$info("- Visualizing overlaps between subcases ...")
                process_overlaps(markers, overlaps, name, each)
            }

        }

        if (!is.null(enriches) && length(enriches) > 0) {
            log$info("- Summarizing enrichments in subcases (by each: {each}) ...")
            if (!is.data.frame(enriches)) {
                each_levels <- names(enriches)
                enriches <- do_call(rbind, lapply(each_levels, function(x) {
                    enrich_df <- enriches[[x]]
                    if (nrow(enrich_df) > 0) {
                        enrich_df[[each]] <- x
                    } else {
                        enrich_df[[each]] <- character(0)  # Empty case
                    }
                    enrich_df
                }))
                enriches[[each]] <- factor(enriches[[each]], levels = each_levels)
            }

            if (length(allenrich_plots) > 0 && !is.null(enriches) && nrow(enriches) > 0) {
                log$info("- Visualizing all enrichments together ...")
                # add other metadata columns if any by mapping groupname
                # only add the metadata columns from object if there is a single value mapped
                metacols <- srtobj@meta.data %>% group_by(!!sym(each)) %>%
                    summarize(across(everything(), ~ n_distinct(.) == 1), .groups = "keep") %>%
                    select(where(~ all(. == TRUE))) %>%
                    colnames()

                if (length(metacols) > 1) {
                    metadf <- srtobj@meta.data[, metacols, drop = FALSE]  %>%
                        distinct(!!sym(each), .keep_all = TRUE)

                    for (col in setdiff(metacols, each)) {
                        if (col %in% colnames(enriches)) {
                            warning("Column name conflict: {col}, adding with suffix '_meta'", immediate. = TRUE)
                            metadf[[paste0(col, "_meta")]] <- metadf[[col]]
                            metadf[[col]] <- NULL
                        }
                    }

                    enriches <- left_join(enriches, metadf, by = each)
                }

                process_allenriches(enriches, allenrich_plots, name, each)
            }
        }

        return(invisible())
    }

    # Let RunSeuratDEAnalysis handle the subset
    case$subset <- subset_
    case$object <- srtobj
    markers <- do_call(RunSeuratDEAnalysis, case)
    case$object <- NULL  # Release memory
    gc()

    subobj <- if (is.null(subset_)) srtobj else filter(srtobj, !!parse_expr(subset_))

    if (is.null(case$ident_1)) {
        all_idents <- unique(as.character(markers[[case$group_by]]))
        enriches <- list()
        for (ident in all_idents) {
            log$info("- {case$group_by}: {ident} ...")
            ident_markers <- markers[markers[[case$group_by]] == ident, , drop = TRUE]
            casename <- paste0(name, "::", paste0(case$group_by, ": ", ident))
            info <- case_info(casename, outdir, create = TRUE)

            attr(ident_markers, "ident_1") <- ident
            enrich <- process_markers(ident_markers, info = info, case = list(
                object = subobj,
                dbs = dbs,
                group_by = case$group_by,
                sigmarkers = sigmarkers,
                enrich_style = enrich_style,
                marker_plots = marker_plots,
                enrich_plots = enrich_plots,
                error = case$error,
                ident = NULL
            ))
            enriches[[ident]] <- enrich
        }

        if (length(allmarker_plots) > 0) {
            log$info("- Visualizing all markers together ...")
            process_allmarkers(
                markers,
                object = subobj,
                comparison_by = case$group_by,
                plotcases = allmarker_plots,
                casename = name,
                groupname = case$group_by,
                subset_by_group = FALSE)
        }

        if (length(overlaps) > 0) {
            log$info("- Visualizing overlaps between subcases ...")
            process_overlaps(markers, overlaps, name, case$group_by)
        }

        if (length(allenrich_plots) > 0) {
            log$info("- Visualizing all enrichments together ...")
            # add other metadata columns if any by mapping groupname
            # only add the metadata columns from object if there is a single value mapped
            metacols <- subobj@meta.data %>% group_by(!!sym(case$group_by)) %>%
                summarize(across(everything(), ~ n_distinct(.) == 1), .groups = "keep") %>%
                select(where(~ all(. == TRUE))) %>%
                colnames()

            if (length(metacols) > 1) {
                metadf <- subobj@meta.data[, metacols, drop = FALSE]  %>%
                    distinct(!!sym(case$group_by), .keep_all = TRUE)

                for (col in setdiff(metacols, case$group_by)) {
                    if (col %in% colnames(enriches[[1]])) {
                        warning("Column name conflict: {col}, adding with suffix '_meta'", immediate. = TRUE)
                        metadf[[paste0(col, "_meta")]] <- metadf[[col]]
                        metadf[[col]] <- NULL
                    }
                }

                for (ne in names(enriches)) {
                    if (!case$group_by %in% colnames(enriches[[ne]])) {
                        enriches[[ne]][[case$group_by]] <- ne
                    }
                    enriches[[ne]] <- left_join(enriches[[ne]], metadf, by = case$group_by)
                }
            }
            enriches <- do_call(rbind, enriches)
            process_allenriches(enriches, allenrich_plots, name, case$group_by)
        }
    } else {
        info <- case_info(name, outdir, create = TRUE)
        enrich <- process_markers(markers, info = info, case = list(
            object = subobj,
            dbs = dbs,
            group_by = case$group_by,
            sigmarkers = sigmarkers,
            enrich_style = enrich_style,
            marker_plots = marker_plots,
            enrich_plots = enrich_plots,
            error = case$error,
            ident = if (is.null(case$ident_2)) case$ident_1 else paste0(case$ident_1, " vs ", case$ident_2)
        ))

        if (!is.null(original_case) && !is.null(cases[[original_case]])) {
            if (nrow(markers) > 0) {
                markers[[each_name]] <- each
            }
            cases[[original_case]]$markers[[each]] <<- markers
            cases[[original_case]]$enriches[[each]] <<- enrich
        }
    }

    invisible()
}

sapply(names(cases), run_case)

reporter$save(joboutdir)
