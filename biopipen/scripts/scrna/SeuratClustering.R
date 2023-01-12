source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(future)
library(tidyr)
library(dplyr)

set.seed(8525)

srtfile = {{in.srtobj | quote}}
rdsfile = {{out.rdsfile | quote}}
envs = {{envs | r}}

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

sobj = readRDS(srtfile)
obj_list = SplitObject(sobj, split.by = "Sample")
rm(sobj)

{% if envs.use_sct -%}
# https://satijalab.org/seurat/articles/integration_rpca.html#performing-integration-on-datasets-normalized-with-sctransform-1
print("- Performing SCTransform on each sample ...")
obj_list <- lapply(X = obj_list, FUN = function(x) {
    print(paste("  Performing SCTransform on sample:", x@meta.data$Sample[1], "..."))
    # # Needed?
    # DefaultAssay(x) <- "RNA"
    args = list_update(envs$SCTransform, list(object = x))
    do_call(SCTransform, args)
})

print("- Running SelectIntegrationFeatures ...")
envs$SelectIntegrationFeatures$object.list = obj_list
features = do_call(SelectIntegrationFeatures, envs$SelectIntegrationFeatures)

print("- Running PrepSCTIntegration ...")
envs$PrepSCTIntegration$object.list = obj_list
envs$PrepSCTIntegration$anchor.features = features
obj_list = do_call(PrepSCTIntegration, envs$PrepSCTIntegration)

print("- Running PCA on each sample ...")
obj_list = lapply(X = obj_list, FUN = RunPCA, features = features)

print("- Running FindIntegrationAnchors ...")
fia_args = list_setdefault(
    envs$FindIntegrationAnchors,
    object.list = obj_list,
    anchor.features = features,
    normalization.method = "SCT",
    reduction = "rpca"
)
anchors = do_call(FindIntegrationAnchors, fia_args)

print("- Running IntegrateData ...")
envs$IntegrateData$anchorset = anchors
id_args = list_setdefault(envs$IntegrateData, normalization.method = "SCT")
obj_list = do_call(IntegrateData, id_args)

{%- else -%}
# https://satijalab.org/seurat/articles/integration_rpca.html
print("- Performing NormalizeData + FindVariableFeatures on each sample ...")
obj_list <- lapply(X = obj_list, FUN = function(x) {
    print(paste("  On sample:", x@meta.data$Sample[1], "..."))
    DefaultAssay(x) <- "RNA"
    args = list_update(envs$NormalizeData, list(object = x))
    x <- do_call(NormalizeData, args)

    args = list_update(envs$FindVariableFeatures, list(object = x))
    do_call(FindVariableFeatures, args)
})


print("- Running SelectIntegrationFeatures ...")
envs$SelectIntegrationFeatures$object.list = obj_list
features = do_call(SelectIntegrationFeatures, envs$SelectIntegrationFeatures)

print("- Running ScaleData + RunPCA on each sample ...")
obj_list <- lapply(X = obj_list, FUN = function(x) {
    print(paste("  On sample:", x@meta.data$Sample[1], "..."))
    args = list_setdefault(envs$ScaleData, object = x, features = features)
    x <- do_call(ScaleData, args)

    RunPCA(x, features = features, verbose = FALSE)
})

print("- Running FindIntegrationAnchors ...")
fia_args = list_setdefault(
    envs$FindIntegrationAnchors,
    object.list = obj_list,
    anchor.features = features,
    reduction = "rpca"
)
anchors = do_call(FindIntegrationAnchors, fia_args)

print("- Running IntegrateData ...")
envs$IntegrateData$anchorset = anchors
obj_list = do_call(IntegrateData, envs$IntegrateData)

DefaultAssay(obj_list) <- "integrated"

envs$ScaleData$object = obj_list
obj_list = do_call(ScaleData, envs$ScaleData)

{%- endif %}

print("- Running RunPCA ...")
envs$RunPCA$object = obj_list
obj_list = do_call(RunPCA, envs$RunPCA)

print("- Running RunUMAP ...")
envs$RunUMAP$object = obj_list
obj_list = do_call(RunUMAP, envs$RunUMAP)

print("- Running FindNeighbors ...")
envs$FindNeighbors$object = obj_list
obj_list = do_call(FindNeighbors, envs$FindNeighbors)

print("- Running FindClusters ...")
envs$FindClusters$object = obj_list
obj_list = do_call(FindClusters, envs$FindClusters)

print("- Saving results ...")
saveRDS(obj_list, file = rdsfile)
