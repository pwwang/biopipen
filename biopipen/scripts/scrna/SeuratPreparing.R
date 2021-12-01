library(Seurat)
library(future)

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

cached_file = file.path(joboutdir, "obj_list.rds")
if (file.exists(cached_file) && file.mtime(cached_file) > file.mtime(metafile)) {
    obj_list = readRDS(cached_file)
    samples = as.character(metadata$Sample)
} else {
    obj_list = list()
    samples = c()
    print("- Reading samples individually ...")
    for (i in seq_len(nrow(metadata))) {
        mdata = metadata[i,, drop=T] # list
        sample = as.character(mdata$Sample)
        print(paste("  Reading", sample, "..."))
        path = as.character(metadata[i, "RNADir", drop=T])
        if (is.na(path) || !is.character(path) || nchar(path) == 0) {
            warning(paste0("No path found for sample: ", sample))
            next
        }

        samples = c(samples, sample)
        exprs = Read10X(data.dir = path)
        obj = CreateSeuratObject(counts=exprs, project=sample)
        obj = SCTransform(object=obj, return.only.var.genes=FALSE, verbose=FALSE)
        obj = RenameCells(obj, add.cell.id = sample)
        obj$percent.mt = PercentageFeatureSet(obj, pattern = "^MT-")
        if (!is.null(envs$qc) && nchar(envs$qc) > 0) {
            obj = subset(obj, subset = {{envs.qc}})
        }
        # Attach meta data
        for (mname in names(mdata)) {
            if (mname %in% c("RNADir", "TCRDir")) { next }
            mdt = mdata[[mname]]
            if (is.factor(mdt)) {
                mdt = levels(mdt)[mdt]
            }
            obj[[mname]] = mdt
        }
        obj_list[[sample]] = obj
    }
    saveRDS(obj_list, cached_file)
}

print('- PrepSCTIntegration ...')
reference.list = unname(obj_list)
features <- SelectIntegrationFeatures(
    object.list = reference.list,
    nfeatures = 3000
)
reference.list <- PrepSCTIntegration(
    object.list = reference.list,
    anchor.features = features
)
reference.list <- lapply(
    X = reference.list,
    FUN = RunPCA,
    features = features
)

print('- FindIntegrationAnchors ...')
# combined.anchors <- FindIntegrationAnchors(reference.list, dims=1:30)
# "long vectors not supported yet"
# global.obj <- IntegrateData(combined.anchors, dims=1:30)

combined.anchors <- FindIntegrationAnchors(
    reference.list,
    normalization.method = "SCT",
    anchor.features = features,
    k.anchor = 20,
    reduction = "rpca",
    reference = 1,
    dims=1:30
)

print('- IntegrateData ...')
global.obj <- IntegrateData(
    combined.anchors,
    normalization.method = "SCT",
    dims=1:30
)
rm(combined.anchors)

print('- RunPCA ...')
global.obj <- RunPCA(global.obj, npcs=30, verbose=FALSE)

print('- RunUMAP ...')
global.obj <- RunUMAP(global.obj, reduction="pca", dims=1:30)

saveRDS(global.obj, rdsfile)
