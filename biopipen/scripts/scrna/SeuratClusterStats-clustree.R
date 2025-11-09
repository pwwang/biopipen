# srtobj, clustrees_defaults, clustrees

log$info("clustrees:")

if (
    (is.null(clustrees) || length(clustrees) == 0) &&
    (is.null(clustrees_defaults$prefix) || isFALSE(clustrees_defaults$prefix))) {
    log$warn("- no case specified, skipping ...")
} else {  # clustrees set or prefix is not empty
    odir = file.path(outdir, "clustrees")
    dir.create(odir, recursive=TRUE, showWarnings=FALSE)

    if ((is.null(clustrees) || length(clustrees) == 0) && isTRUE(clustrees_defaults$prefix)) {
        clustrees <- list()
        for (key in names(srtobj@commands)) {
            if (startsWith(key, "FindClusters") && length(srtobj@commands[[key]]$resolution) > 1) {
                pref <- substring(key, 14)
                if (pref == "") {
                    pref <- biopipen.utils::GetIdentityColumn(srtobj)
                }

                clustrees[[pref]] <- list(prefix = pref)
            }
        }
    }
    if (length(clustrees) == 0) {
        log$warn("- no case found, skipping ...")
    } else {
        reporter$add(
            list(
                kind = "descr",
                content = 'The clustree plots displays clustering results from the Seurat object across different
                resolutions of the clustering algorithm
                (<a target="_blank" href="https://satijalab.org/seurat/reference/findclusters">Seurat::FindClusters</a>).
                Each node represents a cluster, with the resolution levels labeled along the vertical (y) axis.
                The size of each node reflects the number of cells in that cluster. Edges connect clusters between
                adjacent resolutions and indicate how cells transition between clusters as resolution increases.
                The thickness of the edges corresponds to the proportion of shared cells (in_prop) between clusters,
                where darker lines signify a higher overlap (up to 100%). The color of the edges indicates the actual
                number of cells that transitioned between clusters.'
            ),
            h1 = "Clustree plots"
        )

        reports <- list()
        for (name in names(clustrees)) {
            if (is.null(clustrees[[name]]$prefix)) {
                stop(paste0("clustrees: prefix is required for case: ", name))
            }
            log$info("- Case: {name} ...")
            case <- list_update(clustrees_defaults, clustrees[[name]])
            extract_vars(case, "devpars", "more_formats", "save_code")

            prefix <- sub("\\.$", "", case$prefix)
            case$prefix <- paste0(prefix, ".")
            case$object <- srtobj

            command <- srtobj@commands[[paste0("FindClusters.", prefix)]] %||%
                (if(prefix == "seurat_clusters") srtobj@commands$FindClusters else NULL)

            if (is.null(command)) {
                resolution <- substring(colnames(case$x), nchar(case$prefix) + 1)
            } else {
                resolution <- command$resolution
            }
            resolution_used <- resolution[length(resolution)]

            plot_prefix <- file.path(odir, paste0(slugify(prefix), ".clustree"))
            p <- do_call(gglogger::register(ClustreePlot), case)
            save_plot(p, plot_prefix, devpars, formats = c("png", more_formats))

            if (save_code) {
                save_plotcode(p, plot_prefix,
                    setup = c("library(scplotter)", "load('data.RData')", "invisible(list2env(case, envir = .GlobalEnv))"),
                    "case",
                    auto_data_setup = FALSE)
            }
            reports[[length(reports) + 1]] <- reporter$image(
                plot_prefix, more_formats, save_code, kind = "image",
                descr = paste0("Resolutions: ", paste(resolution, collapse = ", "), "; resolution used: ", resolution_used)
            )
        }
        reports$h1 <- "Clustree plots"
        reports$ui <- "table_of_images"
        do_call(reporter$add, reports)
    }
}
