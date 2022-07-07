source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(future)
library(bracer)
library(tidyseurat)

metafile = {{in.metafile | quote}}
rdsfile = {{out.rdsfile | quote}}
joboutdir = {{job.outdir | quote}}
envs = {{envs | r}}

set.seed(8525)
options(future.globals.maxSize = 80000 * 1024^2)
plan(strategy = "multicore", workers = envs$ncores)

metadata = read.table(
    metafile,
    header = TRUE,
    row.names = NULL,
    sep = "\t",
    check.names = FALSE
)

meta_cols = colnames(metadata)
if (!"Sample" %in% meta_cols) {
    stop("Error: Column `Sample` is not found in metafile.")
}
if (!"RNADir" %in% meta_cols) {
    stop("Error: Column `RNADir` is not found in metafile.")
}

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

rename_files = function(e) {
    tmpdatadir = file.path(joboutdir, "renamed", sample)
    if (dir.exists(tmpdatadir)) {
        unlink(tmpdatadir, recursive = TRUE)
    }
    dir.create(tmpdatadir, recursive = TRUE, showWarnings = FALSE)
    barcodefile = Sys.glob(file.path(path, "*barcodes.tsv.gz"))[1]
    file.symlink(
        normalizePath(barcodefile),
        file.path(tmpdatadir, "barcodes.tsv.gz")
    )
    genefile = glob(file.path(path, "*{genes,features}.tsv.gz"))[1]
    file.symlink(
        normalizePath(genefile),
        file.path(tmpdatadir, "features.tsv.gz")
    )
    matrixfile = Sys.glob(file.path(path, "*matrix.mtx.gz"))[1]
    file.symlink(
        normalizePath(matrixfile),
        file.path(tmpdatadir, "matrix.mtx.gz")
    )
    Read10X(data.dir = tmpdatadir)
}

load_sample = function(sample) {
    print(paste("  Loading sample:", sample, "..."))
    mdata = as.data.frame(metadata)[metadata$Sample == sample, , drop=TRUE]
    path = as.character(mdata$RNADir)
    if (is.na(path) || !is.character(path) || nchar(path) == 0) {
        warning(paste0("No path found for sample: ", sample))
        return (NULL)
    }

    # obj_list = list()
    exprs = tryCatch(
        # Read10X requires
        # - barcodes.tsv.gz
        # - genes.tsv.gz
        # - matrix.mtx.gz
        # But sometimes, they are prefixed with sample name
        # e.g.GSM4143656_SAM24345863-ln1.barcodes.tsv.gz
        { Read10X(data.dir = path) },
        error = rename_files
    )
    if ("Gene Expression" %in% names(exprs)) {
        exprs = exprs[["Gene Expression"]]
    }
    obj = CreateSeuratObject(counts=exprs, project=sample)
    # filter the cells that don't have any gene expressions
    cell_exprs = colSums(obj@assays$RNA)
    obj = subset(obj, cells = names(cell_exprs[cell_exprs > 0]))
    # obj = SCTransform(object=obj, return.only.var.genes=FALSE, verbose=FALSE)
    obj = RenameCells(obj, add.cell.id = sample)
    obj$percent.mt = PercentageFeatureSet(obj, pattern = "^MT-")
    if (!is.null(envs$qc) && nchar(envs$qc) > 0) {
        obj = obj |> filter(parse_expr(envs$qc))
    }
    # Attach meta data
    for (mname in names(mdata)) {
        if (mname %in% c("RNADir", "TCRDir")) { next }
        mdt = mdata[[mname]]
        if (is.factor(mdt)) { mdt = levels(mdt)[mdt] }
        obj[[mname]] = mdt
    }
    # obj_list[[sample]] = obj

    # obj_list
    obj
}

# Load data
samples = as.character(metadata$Sample)
cached_file = file.path(joboutdir, "obj_list.rds")
if (file.exists(cached_file) && file.mtime(cached_file) > file.mtime(metafile)) {
    print("- Loading cached data ...")
    obj_list = readRDS(cached_file)
} else {
    print("- Reading samples individually ...")
    obj_list = lapply(samples, load_sample)
    names(obj_list) = samples
    saveRDS(obj_list, file = cached_file)
}

{% if envs.use_sct -%}
# https://satijalab.org/seurat/articles/integration_rpca.html#performing-integration-on-datasets-normalized-with-sctransform-1
print("- Performing SCTransform on each sample ...")
obj_list <- lapply(X = obj_list, FUN = function(x) {
    print(paste("  Performing SCTransform on sample:", x@project.name, "..."))
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
id_args = list_setdefault(envs$IntegrateData, normalization.method = "SCT")
obj_list = do_call(IntegrateData, id_args)

{%- else -%}
# https://satijalab.org/seurat/articles/integration_rpca.html
print("- Performing NormalizeData + FindVariableFeatures on each sample ...")
obj_list <- lapply(X = obj_list, FUN = function(x) {
    print(paste("  On sample:", x@project.name, "..."))
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
    print(paste("  On sample:", x@project.name, "..."))
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

saveRDS(obj_list, rdsfile)
