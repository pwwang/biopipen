source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(future)
library(tidyr)
library(dplyr)
library(digest)

set.seed(8525)

srtfile = {{in.srtobj | quote}}
rdsfile = {{out.rdsfile | quote}}
joboutdir = {{job.outdir | quote}}
envs = {{envs | r: todot="-"}}

options(future.globals.maxSize = 80000 * 1024^2)
plan(strategy = "multicore", workers = envs$ncores)

.expand_dims = function(args, name = "dims") {
    # Expand dims from 30 to 1:30
    if (is.numeric(args[[name]]) && length(args[[name]] == 1)) {
        args[[name]] = 1:args[[name]]
    }
    args
}
envs$FindIntegrationAnchors = .expand_dims(envs$FindIntegrationAnchors)
envs$IntegrateData = .expand_dims(envs$IntegrateData)
envs$RunUMAP = .expand_dims(envs$RunUMAP)
envs$FindNeighbors = .expand_dims(envs$FindNeighbors)

log_info("Reading Seurat object ...")
sobj = readRDS(srtfile)

if (isTRUE(envs$cache)) {
    envs$cache = joboutdir
}

if (is.character(envs$cache) && nchar(envs$cache) > 0) {
    log_info("Obtainning the signature ...")
    envs2 = envs
    envs2$ncores <- NULL
    sig = c(
        capture.output(str(sobj)),
        "\n\n-------------------\n\n",
        capture.output(str(envs2)),
        "\n"
    )
    digested_sig = digest::digest(sig, algo = "md5")
    cached_file = file.path(envs$cache, paste0(digested_sig, ".cached.RDS"))
    if (file.exists(cached_file)) {
        log_info("Using cached results {cached_file}")
        # copy cached file to rdsfile
        file.copy(cached_file, rdsfile, copy.date = TRUE)
        quit()
    } else {
        log_info("Cached results not found, logging the current and cached signatures.")
        log_info("- Current signature:")
        print(sig)
        sigfiles = Sys.glob(file.path(envs$cache, "*.signature.txt"))
        for (sigfile in sigfiles) {
            log_info("- Found cached signature file: {sigfile}")
            cached_sig = readLines(sigfile)
            log_info("- Cached signature:")
            print(cached_sig)
        }
        writeLines(sig, file.path(envs$cache, paste0(digested_sig, ".signature.txt")))
    }
}


obj_list = SplitObject(sobj, split.by = "Sample")
rm(sobj)

# Convert envs$FindIntegrationAnchors$reference to index of given as sample names
samples = unlist(lapply(obj_list, function(x) x@meta.data$Sample[1]))
if (!is.null(envs$FindIntegrationAnchors$reference)) {
    ref = envs$FindIntegrationAnchors$reference
    if (length(ref) == 1) {
        ref = trimws(strsplit(ref, ",")[[1]])
    }
    ref = sapply(ref, function(x) {
        x_int = as.integer(x)
        if (!is.na(x_int)) {
            return(x_int)
        }
        which(samples == x)
    })
    envs$FindIntegrationAnchors$reference = ref
}

{% if envs.use_sct -%}
# ############################
# Using SCT
# https://satijalab.org/seurat/articles/integration_rpca.html#performing-integration-on-datasets-normalized-with-sctransform-1
log_info("########## Using SCT route ##########")
log_info("Performing SCTransform on each sample ...")
obj_list <- lapply(X = obj_list, FUN = function(x) {
    log_info("- On sample: {x@meta.data$Sample[1]} ...")
    # # Needed?
    # DefaultAssay(x) <- "RNA"
    args = list_update(envs$SCTransform, list(object = x))
    do_call(SCTransform, args)
})

log_info("Running SelectIntegrationFeatures ...")
envs$SelectIntegrationFeatures$object.list = obj_list
features = do_call(SelectIntegrationFeatures, envs$SelectIntegrationFeatures)

log_info("Running PrepSCTIntegration ...")
envs$PrepSCTIntegration$object.list = obj_list
envs$PrepSCTIntegration$anchor.features = features
obj_list = do_call(PrepSCTIntegration, envs$PrepSCTIntegration)

log_info("Running PCA on each sample ...")
obj_list = lapply(X = obj_list, FUN = function(x) {
    log_info("- On sample: {x@meta.data$Sample[1]} ...")
    npcs = if (is.null(envs$RunPCA1$npcs)) 50 else envs$RunPCA1$npcs
    args = list_setdefault(
        envs$RunPCA1,
        object = x,
        features = features,
        verbose = FALSE,
        npcs = min(npcs, ncol(x) - 1)
    )
    do_call(RunPCA, args)
})

log_info("Running FindIntegrationAnchors ...")
if (!is.null(envs$FindIntegrationAnchors$reference)) {
    log_info(
        paste(
            "- Using samples as reference:",
            paste(envs$FindIntegrationAnchors$reference, collapse = ", ")
        )
    )
}
fia_args = list_setdefault(
    envs$FindIntegrationAnchors,
    object.list = obj_list,
    anchor.features = features,
    normalization.method = "SCT",
    reduction = "rpca",
    dims = 1:30,
    k.score = 30
)
min_dim = min(unlist(lapply(obj_list, ncol))) - 1
fia_args$dims = 1:min(min_dim, max(fia_args$dims))
fia_args$k.score = min(30, min_dim - 1)
anchors = do_call(FindIntegrationAnchors, fia_args)

log_info("Running IntegrateData ...")
envs$IntegrateData$anchorset = anchors
id_args = list_setdefault(
    envs$IntegrateData,
    normalization.method = "SCT",
    dims = 1:30
)
id_args$dims = 1:min(min_dim, max(id_args$dims))
tryCatch({
    obj_list = do_call(IntegrateData, id_args)
}, error = function(e) {
    msg = ""
    if (grepl("number of items to replace is not a multiple of replacement length", e)) {
        default_kweight = 100
        if (!is.null(envs$IntegrateData$k.weight)) {
            default_kweight = envs$IntegrateData$k.weight
        }
        msg = paste0(
            "It's possible that you have too few cells in some samples, ",
            "causing a small number of anchor cells in the anchorset. \n",
            "  Try changing `k.weight` for `IntegrateData` by setting ",
            "`envs.IntegrateData.k-weight` to a smaller number (it's now ",
            default_kweight, "). \n",
            "  See also https://github.com/satijalab/seurat/issues/6341"
        )
    }
    stop(paste0(msg, "\n", e))
})

{%- else -%}
# ############################
# Using rpca
# https://satijalab.org/seurat/articles/integration_rpca.html
log_info("########## Using rpca route ##########")
log_info("Performing NormalizeData + FindVariableFeatures on each sample ...")
obj_list <- lapply(X = obj_list, FUN = function(x) {
    log_info("- On sample: {x@meta.data$Sample[1]} ...")
    DefaultAssay(x) <- "RNA"
    args = list_update(envs$NormalizeData, list(object = x))
    x <- do_call(NormalizeData, args)

    args = list_update(envs$FindVariableFeatures, list(object = x))
    do_call(FindVariableFeatures, args)
})

log_info("Running SelectIntegrationFeatures ...")
envs$SelectIntegrationFeatures$object.list = obj_list
features = do_call(SelectIntegrationFeatures, envs$SelectIntegrationFeatures)

log_info("Running ScaleData + RunPCA on each sample ...")
obj_list <- lapply(X = obj_list, FUN = function(x) {
    log_info("- On sample: {x@meta.data$Sample[1]} ...")
    args = list_setdefault(envs$ScaleData1, object = x, features = features)
    x <- do_call(ScaleData, args)

    npcs = if (is.null(envs$RunPCA1$npcs)) 50 else envs$RunPCA1$npcs
    args = list_setdefault(
        envs$RunPCA1,
        object = x,
        features = features,
        verbose = FALSE,
        npcs = min(npcs, ncol(x) - 1)
    )
    do_call(RunPCA, args)
})

log_info("Running FindIntegrationAnchors ...")
if (!is.null(envs$FindIntegrationAnchors$reference)) {
    log_info(
        paste(
            "- Using samples as reference:",
            paste(envs$FindIntegrationAnchors$reference, collapse = ", ")
        )
    )
}
fia_args = list_setdefault(
    envs$FindIntegrationAnchors,
    object.list = obj_list,
    anchor.features = features,
    reduction = "rpca",
    dims = 1:30,
    k.score = 30
)
min_dim = min(unlist(lapply(obj_list, ncol))) - 1
fia_args$dims = 1:min(min_dim, max(fia_args$dims))
fia_args$k.score = min(30, min_dim - 1)
anchors = do_call(FindIntegrationAnchors, fia_args)

log_info("Running IntegrateData ...")
envs$IntegrateData$anchorset = anchors
id_args = list_setdefault(envs$IntegrateData, dims = 1:30)
id_args$dims = 1:min(min_dim, max(id_args$dims))
obj_list = do_call(IntegrateData, id_args)

DefaultAssay(obj_list) <- "integrated"

envs$ScaleData$object = obj_list
obj_list = do_call(ScaleData, envs$ScaleData)

{%- endif %}

log_info("Running RunPCA ...")
pca_args = list_setdefault(
    envs$RunPCA,
    object = obj_list,
    npcs = 50
)
pca_args$npcs = min(pca_args$npcs, ncol(obj_list) - 1)
obj_list = do_call(RunPCA, pca_args)

log_info("Running RunUMAP ...")
umap_args = list_setdefault(
    envs$RunUMAP,
    object = obj_list,
    dims = 1:30
)
umap_args$dims = 1:min(max(umap_args$dims), ncol(obj_list) - 1)
obj_list = do_call(RunUMAP, umap_args)

log_info("Running FindNeighbors ...")
envs$FindNeighbors$object = obj_list
obj_list = do_call(FindNeighbors, envs$FindNeighbors)

log_info("Running FindClusters ...")
envs$FindClusters$object = obj_list
obj_list = do_call(FindClusters, envs$FindClusters)

nclusters = length(unique(Idents(obj_list)))
log_info("Identified {nclusters} clusters.")

log_info("Saving results ...")
saveRDS(obj_list, file = rdsfile)

if (is.character(envs$cache) && nchar(envs$cache) > 0) {
    log_info("Caching results ...")
    file.copy(rdsfile, cached_file, overwrite = TRUE)
}
