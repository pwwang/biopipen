library(rlang)
library(dplyr)
library(Seurat)
library(plotthis)
library(biopipen.utils)

log <- get_logger()
reporter <- get_reporter()

srtfile <- {{ in.srtobj | r }}
outdir <- {{ out.outdir | r }}
joboutdir <- {{ job.outdir | r }}

ncores <- {{ envs.ncores | int }}
mutaters <- {{ envs.mutaters | r }}
group.by <- {{ envs["group-by"] | r }}
ident.1 <- {{ envs["ident-1"] | r }}
ident.2 <- {{ envs["ident-2"] | r }}
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
    options(future.globals.maxSize = 80000 * 1024^2)
    plan(strategy = "multicore", workers = ncores)
}

log$info("Reading Seurat object ...")
srtobj <- read_obj(srtfile)
if (!"Identity" %in% colnames(srtobj@meta.data)) {
    srtobj@meta.data$Identity <- Idents(srtobj)
}


if (!is.null(mutaters) && length(mutaters) > 0) {
    log$info("Mutating meta data ...")
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

allmarker_plots <- lapply(allmarker_plots, function(x) {
    list_update(allmarker_plots_defaults, x)
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
    group.by = group.by,
    ident.1 = ident.1,
    ident.2 = ident.2,
    dbs = dbs,
    sigmarkers = sigmarkers,
    enrich_style = enrich_style,
    assay = assay %||% DefaultAssay(srtobj),
    each = each,
    error = error,
    subset = subset,
    allmarker_plots_defaults = allmarker_plots_defaults,
    allmarker_plots = allmarker_plots,
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

    case$group.by <- case$group.by %||% "Identity"

    if (is.null(case$each) || is.na(case$each) || nchar(case$each) == 0 || isFALSE(each)) {
        # single cases, no need to expand
        if (length(case$ident.1) > 0 && length(case$overlaps) > 0) {
            stop("Cannot perform 'overlaps' with a single comparison (ident-1 is set) in case '", name, "'")
        }
        if (length(case$ident.1) > 0 && length(case$allmarker_plots) > 0) {
            stop("Cannot perform 'allmarker_plots' with a single comparison (ident-1 is set) in case '", name, "'")
        }

        case$allmarker_plots <- lapply(
            case$allmarker_plots,
            function(x) { list_update(case$allmarker_plots_defaults, x) }
        )
        case$allmarker_plots_defaults <- NULL

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
        if (length(case$overlaps) > 0 && is.null(case$ident.1)) {
            stop("Cannot perform 'overlaps' analysis with 'each' and without 'ident.1' in case '", name, "'")
        }

        if (length(cases) == 0 && name == "Marker Discovery") {
            name <- case$each
        }

        for (each in eachs) {
            newname <- paste0(name, " - ", each)
            newcase <- case

            newcase$original_case <- name
            newcase$each_name <- case$each
            newcase$each <- each

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
            newcase$overlaps <- NULL
            newcase$overlaps_defaults <- NULL

            outcases[[newname]] <- newcase
        }

        if (length(case$overlaps) > 0 || length(case$allmarker_plots) > 0) {
            ovcase <- case
            ovcase$markers <- list()
            ovcase$allmarker_plots <- lapply(
                ovcase$allmarker_plots,
                function(x) { list_update(ovcase$allmarker_plots_defaults, x) }
            )
            ovcase$allmarker_plots_defaults <- NULL
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

    for (plotname in names(case$marker_plots)) {
        plotargs <- case$marker_plots[[plotname]]
        plotargs$degs <- markers
        rownames(plotargs$degs) <- make.unique(markers$gene)
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

    # Do enrichment analysis
    significant_markers <- unique(sigmarkers$gene)

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

                        attr(p, "height") <- attr(p, "height") / 1.5
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
        })
    }
}

process_allmarkers <- function(markers, plotcases, casename, groupname) {
    name <- paste0(casename, "::", paste0(groupname, " (All Markers)"))
    info <- case_info(name, outdir, create = TRUE)
    for (plotname in names(plotcases)) {
        plotargs <- plotcases[[plotname]]
        plotargs$degs <- markers
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

process_overlaps <- function(markers, ovcases, casename, groupname) {
    name <- paste0(casename, "::", paste0(groupname, ": Overlaps"))
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
    log$info("Case: {name} ...")

    case <- extract_vars(
        case,
        "dbs", "sigmarkers", "allmarker_plots", "marker_plots", "enrich_plots", "overlaps",
        "original_case", "markers", "each_name", "each", "enrich_style",
        allow_nonexisting = TRUE
    )
    if (!is.null(markers)) {  # It is the overlap/allmarker case
        log$info("- Summarizing markers in subcases (by each: {each}) ...")
        # handle the overlaps / allmarkers analysis here
        if (!is.data.frame(markers)) {
            markers <- do_call(rbind, lapply(names(markers), function(x) {
                markers_df <- markers[[x]]
                markers_df[[each]] <- x
                markers_df
            }))
        }
        # gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, diff_pct, <each>

        if (length(allmarker_plots) > 0) {
            log$info("- Visualizing all markers together ...")
            attr(markers, "object") <- srtobj
            attr(markers, "group.by") <- each
            attr(markers, "ident.1") <- NULL
            attr(markers, "ident.2") <- NULL
            process_allmarkers(markers, allmarker_plots, name, each)
        }

        if (length(overlaps) > 0) {
            log$info("- Visualizing overlaps between subcases ...")
            process_overlaps(markers, overlaps, name, each)
        }

        return(invisible())
    }
    case$object <- srtobj
    markers <- do_call(RunSeuratDEAnalysis, case)
    case$object <- NULL
    gc()

    if (is.null(case$ident.1)) {
        all_idents <- unique(as.character(markers[[case$group.by]]))
        for (ident in all_idents) {
            log$info("- {case$group.by}: {ident} ...")
            ident_markers <- markers[markers[[case$group.by]] == ident, , drop = TRUE]
            casename <- paste0(name, "::", paste0(case$group.by, ": ", ident))
            info <- case_info(casename, outdir, create = TRUE)

            attr(ident_markers, "ident.1") <- ident
            process_markers(ident_markers, info = info, case = list(
                dbs = dbs,
                sigmarkers = sigmarkers,
                enrich_style = enrich_style,
                marker_plots = marker_plots,
                enrich_plots = enrich_plots,
                error = case$error,
                ident = NULL
            ))
        }

        if (length(allmarker_plots) > 0) {
            log$info("- Visualizing all markers together ...")
            process_allmarkers(markers, allmarker_plots, name, case$group.by)
        }

        if (length(overlaps) > 0) {
            log$info("- Visualizing overlaps between subcases ...")
            process_overlaps(markers, overlaps, name, case$group.by)
        }
    } else {
        info <- case_info(name, outdir, create = TRUE)
        process_markers(markers, info = info, case = list(
            dbs = dbs,
            sigmarkers = sigmarkers,
            enrich_style = enrich_style,
            marker_plots = marker_plots,
            enrich_plots = enrich_plots,
            error = case$error,
            ident = if (is.null(case$ident.2)) case$ident.1 else paste0(case$ident.1, " vs ", case$ident.2)
        ))

        if (!is.null(original_case)) {
            markers[[each_name]] <- each
            cases[[original_case]]$markers[[each]] <<- markers
        }
    }

    invisible()
}

sapply(names(cases), run_case)

reporter$save(joboutdir)
