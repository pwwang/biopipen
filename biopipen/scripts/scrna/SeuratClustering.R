source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(future)
library(tidyr)
library(dplyr)
library(digest)

set.seed(8525)

srtfile <- {{in.srtobj | quote}}
rdsfile <- {{out.rdsfile | quote}}
joboutdir <- {{job.outdir | quote}}
envs <- {{envs | r: todot="-"}}

if (length(envs$ScaleData) > 0 && length(envs$SCTransform) > 0) {
    stop("Cannot specify both ScaleData and SCTransform")
}

options(str = strOptions(vec.len = 5, digits.d = 5))
options(future.globals.maxSize = 80000 * 1024^2)
plan(strategy = "multicore", workers = envs$ncores)

.expand_dims <- function(args, name = "dims") {
    # Expand dims from 30 to 1:30
    if (is.numeric(args[[name]]) && length(args[[name]] == 1)) {
        args[[name]] <- 1:args[[name]]
    }
    args
}

envs$RunUMAP <- .expand_dims(envs$RunUMAP)
envs$FindNeighbors <- .expand_dims(envs$FindNeighbors)

log_info("Reading Seurat object ...")
sobj <- readRDS(srtfile)

if (isTRUE(envs$cache)) {
    envs$cache <- joboutdir
}

if (is.character(envs$cache) && nchar(envs$cache) > 0) {
    log_info("Obtainning the signature ...")
    envs2 <- envs
    envs2$ncores <- NULL
    sig <- c(
        capture.output(str(sobj)),
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
        log_info("Cached results not found, logging the current and cached signatures.")
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

if (length(envs$ScaleData) > 0) {
    if (DefaultAssay(sobj) == "SCT") {
        stop("SCT assay detected, but ScaleData is specified. Use SCTransform instead.")
    }
    log_info("Running ScaleData ...")
    envs$ScaleData$object <- sobj
    sobj <- do_call(ScaleData, envs$ScaleData)
} else if (length(envs$SCTransform) > 0) {
    if (DefaultAssay(sobj) != "SCT") {
        stop("SCT assay not detected, but SCTransform is specified. Use ScaleData instead.")
    }
    log_info("Running SCTransform ...")
    envs$SCTransform$object <- sobj
    sobj <- do_call(SCTransform, envs$SCTransform)
}

log_info("Running RunUMAP ...")
umap_args <- list_setdefault(
    envs$RunUMAP,
    object = sobj,
    dims = 1:30,
    reduction = sobj@misc$integrated_new_reduction %||% "pca"
)
umap_args$dims <- 1:min(max(umap_args$dims), ncol(sobj) - 1)
sobj <- do_call(RunUMAP, umap_args)

log_info("Running FindNeighbors ...")
envs$FindNeighbors$object <- sobj
if (is.null(envs$FindNeighbors$reduction)) {
    envs$FindNeighbors$reduction <- sobj@misc$integrated_new_reduction %||% "pca"
}
sobj <- do_call(FindNeighbors, envs$FindNeighbors)

log_info("Running FindClusters ...")
if (is.null(envs$FindClusters$random.seed)) {
    envs$FindClusters$random.seed <- 8525
}
resolution <- envs$FindClusters$resolution
if (is.character(resolution)) {
    if (grepl(",", resolution)) {
        resolution <- as.numeric(trimws(unlist(strsplit(resolution, ","))))
    } else {
        resolution <- as.numeric(resolution)
    }
}
if (is.null(resolution) || length(resolution) == 1) {
    envs$FindClusters$resolution <- resolution
    envs$FindClusters$object <- sobj
    sobj <- do_call(FindClusters, envs$FindClusters)
    levels(sobj$seurat_clusters) <- paste0("c", as.numeric(levels(sobj$seurat_clusters)) + 1)
    Idents(sobj) <- "seurat_clusters"
    ident_table <- table(sobj$seurat_clusters)
    log_info("- Found {length(ident_table)} clusters:")
    print(ident_table)
} else {
    log_info("- Multiple resolutions detected ...")
    res_key <- NULL
    for (res in resolution) {
        findclusters_args <- envs$FindClusters
        findclusters_args$resolution <- res
        findclusters_args$object <- sobj
        sobj <- do_call(FindClusters, findclusters_args)
        levels(sobj$seurat_clusters) <- paste0("c", as.numeric(levels(sobj$seurat_clusters)) + 1)
        res_key <- paste0("seurat_clusters_", res)
        sobj[[res_key]] <- sobj$seurat_clusters
        ident_table <- table(sobj[[res_key]])
        log_info("- Found {length(ident_table)} at resolution: {res}:")
        print(ident_table)
    }
}

if (DefaultAssay(sobj) == "SCT") {
    # https://github.com/satijalab/seurat/issues/6968
    log_info("Running PrepSCTFindMarkers ...")
    sobj <- PrepSCTFindMarkers(sobj)
}

log_info("Saving results ...")
saveRDS(sobj, file = rdsfile)

if (is.character(envs$cache) && nchar(envs$cache) > 0) {
    log_info("Caching results ...")
    file.copy(rdsfile, cached_file, overwrite = TRUE)
}
