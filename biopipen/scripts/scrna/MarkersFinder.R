library(rlang)
library(dplyr)
library(Seurat)
library(plotthis)
library(biopipen.utils)

log <- get_logger()
reporter <- get_reporter()

srtfile <- {{ in.srtobj | quote }}
outdir <- {{ out.outdir | quote }}
joboutdir <- {{ job.outdir | quote }}
ncores <- {{ envs.ncores | int }}
mutaters <- {{ envs.mutaters | r }}
ident.1 <- {{ envs["ident-1"] | r }}
ident.2 <- {{ envs["ident-2"] | r }}
group.by <- {{ envs["group-by"] | r }}
each <- {{ envs.each | r }}
prefix_each <- {{ envs.prefix_each | r }}
prefix_group <- {{ envs.prefix_group | r }}
assay <- {{ envs.assay | r }}
subset <- {{ envs.subset | r }}
error <- {{ envs.error | r }}
site <- {{ envs.site | r }}
rest <- {{ envs.rest | r: todot="-" }}
dbs <- {{ envs.dbs | r }}
sigmarkers <- {{ envs.sigmarkers | r }}
cache <- {{ envs.cache | r }}
allmarker_plots_defaults <- {{ envs.allmarker_plots_defaults | r }}
allmarker_plots <- {{ envs.allmarker_plots | r }}
marker_plots_defaults <- {{ envs.marker_plots_defaults | r }}
marker_plots <- {{ envs.marker_plots | r }}
enrich_plots_defaults <- {{ envs.enrich_plots_defaults | r }}
enrich_plots <- {{ envs.enrich_plots | r }}
cases <- {{ envs.cases | r: todot="-", skip=1 }}
overlaps_defaults <- {{ envs.overlaps_defaults | r }}
overlaps <- {{ envs.overlaps | r }}

if (isTRUE(cache)) { cache <- joboutdir }

set.seed(8525)
if (ncores > 1) {
    options(future.globals.maxSize = 80000 * 1024^2)
    plan(strategy = "multicore", workers = ncores)
}

log$info("Reading Seurat object ...")
srtobj <- readRDS(srtfile)

if (!is.null(mutaters) && length(mutaters) > 0) {
    log$info("Mutating meta data ...")
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

defaults <- list(
    ident.1 = ident.1,
    ident.2 = ident.2,
    group.by = group.by,
    each = each,
    prefix_each = prefix_each,
    prefix_group = prefix_group,
    dbs = dbs,
    assay = assay %||% DefaultAssay(srtobj),
    subset = subset,
    error = error,
    site = site,
    sigmarkers = sigmarkers,
    allmarker_plots = allmarker_plots,
    marker_plots = marker_plots,
    enrich_plots = enrich_plots,
    cache = cache,
    rest = rest
)

log$info("Expanding cases ...")

post_casing <- function(name, case) {
    outcases <- list()
    no_each <- is.null(case$each) || is.na(case$each) || nchar(case$each) == 0

    if (no_each) {
        # single cases, no need to expand
        case$allmarker_plots <- lapply(
            case$allmarker_plots,
            function(x) { list_update(allmarker_plots_defaults, x) }
        )
        case$marker_plots <- lapply(
            case$marker_plots,
            function(x) { list_update(marker_plots_defaults, x) }
        )
        case$enrich_plots <- lapply(
            case$enrich_plots,
            function(x) { list_update(enrich_plots_defaults, x) }
        )
        outcases[[name]] <- case
    } else {  # !no_each
        if (!is.null(case$subset)) {
            sobj <- srtobj %>% filter(!!parse_expr(case$subset))
        } else {
            sobj <- srtobj
        }

        eachs <- sobj@meta.data %>% pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
        case_1 <- case
        for (each in eachs) {
            each_name <- ifelse(case_1$prefix_each, paste0(case_1$each, " - ", each), each)
            if (!is.null(case_1$ident.1)) {
                # Make name a section
                key <- paste0(name, "::", each_name)
            } else {
                key <- paste0(name, ": ", each_name)
            }
            if (!is.null(case$subset)) {
                case_1$subset <- paste0(case$subset, " & `", case_1$each, "` == '", each, "'")
            } else {
                case_1$subset <- paste0("`", case_1$each, "` == '", each, "'")
            }
            case_1$allmarker_plots <- lapply(
                case_1$allmarker_plots,
                function(x) { list_update(allmarker_plots_defaults, x) }
            )
            case_1$marker_plots <- lapply(
                case_1$marker_plots,
                function(x) { list_update(marker_plots_defaults, x) }
            )
            case_1$enrich_plots <- lapply(
                case_1$enrich_plots,
                function(x) { list_update(enrich_plots_defaults, x) }
            )
            outcases[[key]] <- case_1
        }
    }
    outcases
}
cases <- expand_cases(cases, defaults, post_casing)

# Checking the overlapping cases
case_markers <- list()
if (length(overlaps) > 0) {
    log$info("Checking overlapping cases ...")
    overlaps <- expand_cases(overlaps, overlaps_defaults)
    for (ovname in names(overlaps)) {
        ov <- overlaps[[ovname]]
        # check the existence of the cases
        for (case in ov$cases) {
            if (is.null(cases[[case]])) {
                stop(paste0("Case '", case, "' not found in the cases for overlapping case '", ovname, "'"))
            }
        }
        if (length(ov$cases) < 2) {
            stop("Overlapping cases must have at least 2 cases for overlapping case '", ovname, "'")
        }
        for (case in ov$cases) {
            case_markers[[case]] <- TRUE
        }
        if (identical(ov$venn$enabled, "auto")) {
            overlaps[[ovname]]$venn$enabled <- length(ov$cases) <= 5
        }
    }
}

log$info("Running cases ...")

process_markers <- function(markers, info, case) {
    # Save markers
    write.table(markers, file.path(info$prefix, "markers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    reporter$add2(
        list(
            name = "Table",
            contents = list(list(kind = "table", src = file.path(info$prefix, "markers.tsv"), data = list(nrows = 100)))
        ),
        hs = c(info$section, info$name),
        hs2 = "Markers",
        ui = "tabs"
    )

    for (plotname in names(case$marker_plots)) {
        plotargs <- case$marker_plots[[plotname]]
        plotargs$degs <- markers
        plotargs$outprefix <- file.path(info$prefix, paste0("markers.", slugify(plotname)))
        do_call(VizDEGs, plotargs)
        reporter$add2(
            list(
                name = plotname,
                contents = list(reporter$image(plotargs$outprefix, plotargs$more_formats, plotargs$save_code))),
            hs = c(info$section, info$name),
            hs2 = "Markers",
            ui = "tabs"
        )
    }

    # Do enrichment analysis
    tryCatch({
        enrich <- RunEnrichment(
            markers, deg = case$sigmarkers, dbs = case$dbs, cache = case$cache,
            error = TRUE, site = case$site)

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
                    plotargs$enrich <- enrich[enrich$db == db, , drop = FALSE]
                    plotargs$outprefix <- file.path(info$prefix, paste0("enrich.", slugify(db), ".", slugify(plotname)))

                    do_call(VizEnrich, plotargs)

                    plots[[length(plots) + 1]] <- reporter$image(plotargs$outprefix, plotargs$more_formats, plotargs$save_code)
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

run_case <- function(name) {
    case <- cases[[name]]
    log$info("- Case: {name} ...")

    args <- case$rest %||% list()
    args$object <- srtobj
    args$group.by <- case$group.by
    args$ident.1 <- case$ident.1
    args$ident.2 <- case$ident.2
    args$cache <- case$cache
    args$assay <- case$assay
    args$error <- case$error
    args$subset <- case$subset

    markers <- do_call(RunSeuratDEAnalysis, args)
    if (isTRUE(case_markers[[name]])) {
        case_markers[[name]] <<- markers
    }
    if (is.null(case$ident.1)) {
        if (!is.null(case_markers[[name]])) {
            stop("Case '", name, "' for overlapping analysis must have 'ident.1' defined")
        }
        all_idents <- unique(markers[[case$group.by]])
        # Visualize all markers
        if (length(case$allmarker_plots) > 0) {
            log$info("  Visualizing all markers ...")
            casename <- paste0(name, "::", ifelse(case$prefix_group, paste0(case$group.by, " - All Markers"), "All Markers"))
            info <- case_info(casename, outdir, create = TRUE)
            for (plotname in names(case$allmarker_plots)) {
                plotargs <- case$allmarker_plots[[plotname]]
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
        for (ident in all_idents) {
            log$info("  {case$group.by}: {ident} ...")
            ident_markers <- markers[markers[[case$group.by]] == ident, , drop = TRUE]
            casename <- paste0(name, "::", ifelse(case$prefix_group, paste0(case$group.by, " - ", ident), ident))
            info <- case_info(casename, outdir, create = TRUE)

            process_markers(ident_markers, info = info, case = case)
        }
    } else {
        info <- case_info(name, outdir, create = TRUE)
        process_markers(markers, info = info, case = case)
    }
}

sapply(names(cases), run_case)

if (length(overlaps) > 0) {
    log$info("Running overlapping cases ...")

    run_overlap <- function(ovname) {
        ov <- overlaps[[ovname]]
        ov$sigmarkers <- ov$sigmarkers %||% sigmarkers
        log$info("- Overlapping case: {ovname} ...")
        markers <- lapply(ov$cases, function(case) {
            case_markers[[case]] %>% filter(!!parse_expr(ov$sigmarkers)) %>%
                pull("gene") %>% unique()
        })
        names(markers) <- ov$cases
        info <- case_info(paste0("OVERLAPPING::", ovname), outdir, create = TRUE)

        if (ov$venn$enabled) {
            venn <- extract_vars(ov$venn, "enabled", "more_formats", "save_code", "devpars")
            venn$data <- markers
            venn$in_form <- "list"
            prefix <- file.path(info$prefix, "venn")
            p <- do_call(gglogger::register(VennDiagram), venn)
            save_plot(p, prefix, devpars, formats = c("png", more_formats))
            if (save_code) {
                save_plotcode(
                    p, prefix,
                    c("library(plotthis)", "load('data.RData')", "invisible(list2env(venn, .GlobalEnv))"),
                    "venn",
                    auto_data_setup = FALSE)
            }

            reporter$add2(
                list(
                    name = "Venn Diagram",
                    contents = list(reporter$image(prefix, more_formats, save_code))
                ),
                hs = c(info$section, info$name),
                ui = "tabs"
            )
        }

        if (ov$upset$enabled) {
            upset <- extract_vars(ov$upset, "enabled", "more_formats", "save_code", "devpars")
            upset$data <- markers
            upset$in_form <- "list"
            prefix <- file.path(info$prefix, "upset")
            p <- do_call(gglogger::register(UpsetPlot), upset)
            save_plot(p, prefix, devpars, formats = c("png", more_formats))
            if (save_code) {
                save_plotcode(
                    p, prefix,
                    c("library(plotthis)", "load('data.RData')", "invisible(list2env(upset, .GlobalEnv))"),
                    "upset",
                    auto_data_setup = FALSE)
            }

            reporter$add2(
                list(
                    name = "UpSet Plot",
                    contents = list(reporter$image(prefix, more_formats, save_code))
                ),
                hs = c(info$section, info$name),
                ui = "tabs"
            )
        }

    }

    sapply(names(overlaps), run_overlap)
}

reporter$save(joboutdir)
