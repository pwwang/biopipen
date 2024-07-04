# srtobj, clustrees_defaults, clustrees
log_info("clustrees:")
if (
    (is.null(clustrees) || length(clustrees) == 0) &&
    (is.null(clustrees_defaults$prefix) || clustrees_defaults$prefix == "")) {
    log_warn("- no cases, skipping intentionally ...")
} else {  # clustrees set or prefix is not empty
    library(clustree)
    odir = file.path(outdir, "clustrees")
    dir.create(odir, recursive=TRUE, showWarnings=FALSE)

    if ((is.null(clustrees) || length(clustrees) == 0) && clustrees_defaults$prefix == "_auto") {
        clustrees <- list()
        for (key in names(srtobj@commands)) {
            if (startsWith(key, "FindClusters") && length(srtobj@commands[[key]]$resolution) > 1) {
                pref <- substring(key, 14)
                if (pref == "") {
                    pref <- "seurat_clusters"
                }

                clustrees[[pref]] <- list(prefix = pref)
            }
        }
    }
    if (length(clustrees) == 0) {
        log_warn("- no cases found, skipping ...")
    } else {
        reports <- list()
        for (name in names(clustrees)) {
            if (is.null(clustrees[[name]]$prefix)) {
                stop(paste0("clustrees: prefix is required for case: ", name))
            }
            case <- list_update(clustrees_defaults, clustrees[[name]])

            devpars <- case$devpars
            devpars$width <- devpars$width %||% clustrees_defaults$devpars$width %||% 800
            devpars$height <- devpars$height %||% clustrees_defaults$devpars$height %||% 1000
            devpars$res <- devpars$res %||% clustrees_defaults$devpars$res %||% 100
            case$devpars <- NULL
            prefix <- sub("\\.$", "", case$prefix)
            log_info("- Case: {name} ...")
            case$prefix <- paste0(prefix, ".")
            case$x <- srtobj@meta.data %>% select(starts_with(case$prefix))
            case$x <- case$x[complete.cases(case$x), , drop = FALSE]

            command <- srtobj@commands[[paste0("FindClusters.", prefix)]] %||%
                (if(prefix == "seurat_clusters") srtobj@commands$FindClusters else NULL)

            clustree_file <- file.path(odir, paste0(prefix, ".clustree.png"))
            png(clustree_file, width = devpars$width, height = devpars$height, res = devpars$res)
            p <- do_call(clustree, case)
            print(p)
            dev.off()

            if (is.null(command)) {
                resolution <- substring(colnames(case$x), nchar(case$prefix) + 1)
            } else {
                resolution <- command$resolution
            }
            resolution_used <- resolution[length(resolution)]

            reports[[length(reports) + 1]] <- list(
                kind = "table_image",
                src = clustree_file,
                name = name,
                descr = paste0("Resolutions: ", paste(resolution, collapse = ", "), "; resolution used: ", resolution_used)
            )
        }
        reports$h1 <- "Clustree plots"
        reports$ui <- "table_of_images"
        do.call(add_report, reports)
    }
}
