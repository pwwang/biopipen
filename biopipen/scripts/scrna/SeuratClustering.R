source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/caching.R")

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

if (isTRUE(envs$cache)) { envs$cache <- joboutdir }
if (length(envs$cache) > 1) {
    log_warn("Multiple cache directories (envs.cache) detected, using the first one.")
    envs$cache <- envs$cache[1]
}
sobj_sig <- capture.output(str(sobj))
dig_sig <- digest::digest(sobj_sig, algo = "md5")
dig_sig <- substr(dig_sig, 1, 8)
cache_dir <- NULL
if (is.character(envs$cache)) {
    cache_dir <- file.path(envs$cache, paste0(dig_sig, ".seurat_cache"))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(sobj_sig, file.path(cache_dir, "signature.txt"))
}

if (length(envs$ScaleData) > 0) {
    if (DefaultAssay(sobj) == "SCT") {
        stop("SCT assay detected, but ScaleData is specified. Use SCTransform instead.")
    }
    cached <- get_cached(envs$ScaleData, "ScaleData", cache_dir)
    if (is.null(cached$data)) {
        log_info("Running ScaleData ...")
        envs$ScaleData$object <- sobj
        sobj <- do_call(ScaleData, envs$ScaleData)
        cached$data <- list(assay = sobj@assays$RNA, commands = sobj@commands)
        save_to_cache(cached, "ScaleData", cache_dir)
    } else {
        log_info("Loading cached ScaleData ...")
        sobj@assays$RNA <- cached$data$assay
        sobj@commands <- cached$data$commands
        DefaultAssay(sobj) <- "RNA"
    }
} else if (length(envs$SCTransform) > 0) {
    if (DefaultAssay(sobj) != "SCT") {
        stop("SCT assay not detected, but SCTransform is specified. Use ScaleData instead.")
    }
    cached <- get_cached(envs$SCTransform, "SCTransform", cache_dir)
    asssay <- envs$SCTransform$new.assay.name %||% "SCT"
    if (is.null(cached$data)) {
        log_info("Running SCTransform ...")
        envs$SCTransform$object <- sobj
        sobj <- do_call(SCTransform, envs$SCTransform)
        cached$data <- list(assay = sobj@assays$SCT, commands = sobj@commands)
        save_to_cache(cached, "SCTransform", cache_dir)
    } else {
        log_info("Loading cached SCTransform ...")
        sobj@assays[[assay]] <- cached$data$assay
        sobj@commands <- cached$data$commands
        DefaultAssay(sobj) <- assay
    }
}

cached <- get_cached(envs$RunUMAP, "RunUMAP", cache_dir)
reduc_name <- envs$RunUMAP$reduction.name %||% "umap"
if (is.null(cached$data)) {
    log_info("Running RunUMAP ...")
    umap_args <- list_setdefault(
        envs$RunUMAP,
        object = sobj,
        dims = 1:30,
        reduction = sobj@misc$integrated_new_reduction %||% "pca"
    )
    ncells <- ncol(sobj)
    umap_args$dims <- 1:min(max(umap_args$dims), ncells - 1)
    umap_method <- envs$RunUMAP$umap.method %||% "uwot"
    if (umap_method == "uwot" && is.null(envs$RunUMAP$n.neighbors)) {
        # https://github.com/satijalab/seurat/issues/4312
        umap_args$n.neighbors <- min(ncells - 1, 30)
    }
    sobj <- do_call(RunUMAP, umap_args)
    cached$data <- list(reduc = sobj@reductions[[reduc_name]], commands = sobj@commands)
    save_to_cache(cached, "RunUMAP", cache_dir)
} else {
    log_info("Loading cached RunUMAP ...")
    sobj@reductions[[reduc_name]] <- cached$data$reduc
    sobj@commands <- cached$data$commands
}

cached <- get_cached(envs$FindNeighbors, "FindNeighbors", cache_dir)
if (is.null(cached$data)) {
    log_info("Running FindNeighbors ...")
    envs$FindNeighbors$object <- sobj
    envs$FindNeighbors$reduction <- sobj@misc$integrated_new_reduction %||% "pca"
    sobj <- do_call(FindNeighbors, envs$FindNeighbors)
    cached$data <- list(graphs = sobj@graphs, commands = sobj@commands)
    save_to_cache(cached, "FindNeighbors", cache_dir)
} else {
    log_info("Loading cached FindNeighbors ...")
    sobj@graphs <- cached$data$graphs
    sobj@commands <- cached$data$commands
}

envs$FindClusters$random.seed <- envs$FindClusters$random.seed %||% 8525
resolution <- envs$FindClusters$resolution %||% 0.8
if (is.character(resolution)) {
    if (grepl(",", resolution)) {
        resolution <- as.numeric(trimws(unlist(strsplit(resolution, ","))))
    } else {
        resolution <- as.numeric(resolution)
    }
}

for (res in resolution) {
    envs$FindClusters$resolution <- res
    cached <- get_cached(envs$FindClusters, paste0("FindClusters_", res), cache_dir)
    res_key <- paste0("seurat_clusters_", res)
    if (is.null(cached$data)) {
        log_info("Running FindClusters at resolution: {res} ...")
        envs$FindClusters$object <- sobj
        sobj <- do_call(FindClusters, envs$FindClusters)
        levels(sobj$seurat_clusters) <- paste0("c", as.numeric(levels(sobj$seurat_clusters)) + 1)
        sobj[[res_key]] <- sobj$seurat_clusters
        Idents(sobj) <- "seurat_clusters"
        cached$data <- list(clusters = sobj$seurat_clusters, commands = sobj@commands)
        save_to_cache(cached, paste0("FindClusters_", res), cache_dir)
    } else {
        log_info("Loading cached FindClusters at resolution: {res} ...")
        sobj@commands <- cached$data$commands
        sobj[[res_key]] <- cached$data$clusters
        sobj$seurat_clusters <- cached$data$clusters
        Idents(sobj) <- "seurat_clusters"
    }
    ident_table <- table(Idents(sobj))
    log_info("- Found {length(ident_table)} clusters")
    print(ident_table)
    cat("\n")
}

if (DefaultAssay(sobj) == "SCT") {
        # https://github.com/satijalab/seurat/issues/6968
    log_info("Running PrepSCTFindMarkers ...")
    sobj <- PrepSCTFindMarkers(sobj)
}

log_info("Saving results ...")
saveRDS(sobj, file = rdsfile)
