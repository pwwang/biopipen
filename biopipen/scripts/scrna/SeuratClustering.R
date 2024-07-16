{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "caching.R" | source_r }}

library(Seurat)
library(future)
library(rlang)
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

{{ biopipen_dir | joinpaths: "scripts", "scrna", "SeuratClustering-common.R" | source_r }}

envs$RunUMAP <- expand_dims(envs$RunUMAP)
envs$FindNeighbors <- expand_dims(envs$FindNeighbors)

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

sobj <- run_transformation(sobj)
sobj <- run_umap(sobj)
sobj <- run_findneighbors(sobj)
sobj <- run_findclusters(sobj)
sobj <- run_prepsctfindmarkers(sobj)

log_info("Saving results ...")
saveRDS(sobj, file = rdsfile)
