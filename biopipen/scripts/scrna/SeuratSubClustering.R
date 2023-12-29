source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(future)
library(rlang)
library(tidyr)
library(dplyr)
library(tidyseurat)
library(digest)

set.seed(8525)

srtfile <- {{in.srtobj | quote}}
rdsfile <- {{out.rdsfile | quote}}
joboutdir <- {{job.outdir | quote}}
envs <- {{envs | r: todot = "-"}}

options(str = strOptions(vec.len = 5, digits.d = 5))
options(future.globals.maxSize = 80000 * 1024^2)
plan(strategy = "multicore", workers = envs$ncores)

.expand_dims <- function(args, name = "dims") {
    # Expand dims from 30 to 1:30
    if (!is.null(args) && is.numeric(args[[name]]) && length(args[[name]] == 1)) {
        args[[name]] <- 1:args[[name]]
    }
    args
}

envs$RunUMAP <- .expand_dims(envs$RunUMAP)
envs$FindNeighbors <- .expand_dims(envs$FindNeighbors)

log_info("Reading Seurat object ...")
srtobj <- readRDS(srtfile)

if (isTRUE(envs$cache)) {
    envs$cache <- joboutdir
}

if (is.character(envs$cache) && nchar(envs$cache) > 0) {
    log_info("Obtainning the signature ...")
    envs2 <- envs
    envs2$ncores <- NULL
    sig <- c(
        capture.output(str(srtobj)),
        "\n\n-------------------\n\n",
        capture.output(str(envs2)),
        "\n"
    )
    digested_sig <- digest::digest(sig, algo = "md5")
    cached_file <- file.path(envs$cache, paste0(digested_sig, ".cached.RDS"))
    if (file.exists(cached_file)) {
        log_info("Using cached results {cached_file}")
        # copy cached file to rdsfile
        file.copy(cached_file, rdsfile, copy.date = TRUE)
        quit()
    } else {
        log_info("Cached results not found.")
        log_info("- Current signature: {digested_sig}")
        # print(sig)
        # sigfiles <- Sys.glob(file.path(envs$cache, "*.signature.txt"))
        # for (sigfile in sigfiles) {
        #     log_info("- Found cached signature file: {sigfile}")
        #     cached_sig <- readLines(sigfile)
        #     log_info("- Cached signature:")
        #     print(cached_sig)
        # }
        writeLines(sig, file.path(envs$cache, paste0(digested_sig, ".signature.txt")))
    }
}

if (!is.null(envs$mutaters) && length(envs$mutaters) > 0) {
    log_info("Mutating Seurat object ...")
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

if (length(envs$cases) == 0) {
    envs$cases <- list(subcluster = list())
}

for (key in names(envs$cases)) {
    log_info("")
    log_info("Running case '{key}' ...")
    log_info("===========================================")
    case <- envs$cases[[key]]
    case$RunUMAP <- .expand_dims(case$RunUMAP)
    case$FindNeighbors <- .expand_dims(case$FindNeighbors)

    case <- list_update(
        list(
            subset = envs$subset,
            RunUMAP = envs$RunUMAP,
            FindNeighbors = envs$FindNeighbors,
            FindClusters = envs$FindClusters
        ),
        case
    )

    if (is.null(case$subset) || length(case$subset) == 0) {
        stop(paste0("`subset` for case '", key, "' is empty."))
    }

    log_info("- Subsetting ...")
    sobj <- srtobj %>% filter(!!parse_expr(case$subset))

    log_info("- Running RunUMAP ...")
    umap_args <- list_setdefault(
        case$RunUMAP,
        object = sobj,
        dims = 1:30,
        reduction = sobj@misc$integrated_new_reduction %||% "pca"
    )
    umap_args$dims <- 1:min(max(umap_args$dims), ncol(sobj) - 1)
    sobj <- do_call(RunUMAP, umap_args)

    log_info("- Running FindNeighbors ...")
    case$FindNeighbors$object <- sobj
    if (is.null(case$FindNeighbors$reduction)) {
        case$FindNeighbors$reduction <- sobj@misc$integrated_new_reduction %||% "pca"
    }
    sobj <- do_call(FindNeighbors, case$FindNeighbors)

    log_info("- Running FindClusters ...")
    if (is.null(case$FindClusters$random.seed)) {
        case$FindClusters$random.seed <- 8525
    }
    resolution <- case$FindClusters$resolution
    if (is.character(resolution)) {
        if (grepl(",", resolution)) {
            resolution <- as.numeric(trimws(unlist(strsplit(resolution, ","))))
        } else {
            resolution <- as.numeric(resolution)
        }
    }
    if (is.null(resolution) || length(resolution) == 1) {
        case$FindClusters$resolution <- resolution
        case$FindClusters$object <- sobj
        sobj <- do_call(FindClusters, case$FindClusters)
        levels(sobj$seurat_clusters) <- paste0("s", as.numeric(levels(sobj$seurat_clusters)) + 1)
        Idents(sobj) <- "seurat_clusters"
        sobj[[key]] <- sobj$seurat_clusters
        ident_table <- table(sobj[[key]])
        log_info("- Found {length(ident_table)} clusters:")
        print(ident_table)
        cat("\n")

        log_info("- Updating meta.data with subclusters...")
        srtobj <- AddMetaData(srtobj, metadata = sobj@meta.data[, key, drop = FALSE])
        srtobj[[paste0("sub_umap_", key)]] <- sobj@reductions$umap
    } else {
        log_info("- Multiple resolutions detected ...")
        log_info("")
        metadata <- NULL
        for (res in resolution) {
            findclusters_args <- case$FindClusters
            findclusters_args$resolution <- res
            findclusters_args$object <- sobj
            sobj1 <- do_call(FindClusters, findclusters_args)
            res_key <- paste0(key, "_", res)
            levels(sobj1$seurat_clusters) <- paste0("s", as.numeric(levels(sobj1$seurat_clusters)) + 1)
            Idents(sobj1) <- "seurat_clusters"
            sobj1[[res_key]] <- sobj1$seurat_clusters
            ident_table <- table(sobj1[[res_key]])
            log_info("- Found {length(ident_table)} at resolution: {res}:")
            print(ident_table)
            cat("\n")

            log_info("- Updating meta.data with subclusters...")
            metadata <- sobj1@meta.data[, res_key, drop = FALSE]
            srtobj <- AddMetaData(srtobj, metadata = metadata)
            srtobj[[paste0("sub_umap_", res_key)]] <- sobj1@reductions$umap
        }
        srtobj <- AddMetaData(srtobj, metadata = metadata, col.name = key)
        srtobj[[paste0("sub_umap_", key)]] <- sobj1@reductions$umap
    }
}

log_info("Saving results ...")
saveRDS(srtobj, file = rdsfile)

if (is.character(envs$cache) && nchar(envs$cache) > 0) {
    log_info("Caching results to {cached_file} ...")
    invisible(file.copy(rdsfile, cached_file, overwrite = TRUE))
}
