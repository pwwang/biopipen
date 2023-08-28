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
if (!"RNAData" %in% meta_cols) {
    stop("Error: Column `RNAData` is not found in metafile.")
}


rename_files = function(e, sample, path) {
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
    path = as.character(mdata$RNAData)
    if (is.na(path) || !is.character(path) || nchar(path) == 0) {
        warning(paste0("No path found for sample: ", sample))
        return (NULL)
    }

    # obj_list = list()
    if (dir.exists(path)) {
        exprs = tryCatch(
            # Read10X requires
            # - barcodes.tsv.gz
            # - genes.tsv.gz
            # - matrix.mtx.gz
            # But sometimes, they are prefixed with sample name
            # e.g.GSM4143656_SAM24345863-ln1.barcodes.tsv.gz
            { Read10X(data.dir = path) },
            error = function(e) rename_files(e, sample, path)
        )
    } else {
        exprs = Read10X_h5(path)
    }
    if ("Gene Expression" %in% names(exprs)) {
        exprs = exprs[["Gene Expression"]]
    }
    obj = CreateSeuratObject(counts=exprs, project=sample)
    # filter the cells that don't have any gene expressions
    cell_exprs = colSums(obj@assays$RNA)
    obj = subset(obj, cells = names(cell_exprs[cell_exprs > 0]))
    # obj = SCTransform(object=obj, return.only.var.genes=FALSE, verbose=FALSE)
    obj = RenameCells(obj, add.cell.id = sample)
    # Attach meta data
    for (mname in names(mdata)) {
        if (mname %in% c("RNAData", "TCRData")) { next }
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

print("- Reading samples individually ...")
obj_list = lapply(samples, load_sample)

print("- Merging samples ...")
if (length(obj_list) >= 2) {
    y = c()
    for (i in 2:length(obj_list)) y = c(y, obj_list[[i]])
    sobj = merge(obj_list[[1]], y)
} else {
    sobj = obj_list[[1]]
}

print("- Adding metadata for QC ...")
sobj$percent.mt = PercentageFeatureSet(sobj, pattern = "^MT-")
sobj$percent.ribo = PercentageFeatureSet(sobj, pattern = "^RP[SL]")
sobj$percent.hb = PercentageFeatureSet(sobj, pattern = "^HB[^(P)]")
sobj$percent.plat = PercentageFeatureSet(sobj, pattern = "PECAM1|PF4")

print("- Plotting QC metrics before filtering ...")
feats = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat")
bqcdir = file.path(joboutdir, "before-qc")
dir.create(bqcdir, recursive = TRUE, showWarnings = FALSE)
for (feat in feats) {
    png(file.path(bqcdir, paste0(feat, ".png")), width = 800 + length(samples) * 15, height = 600, res = 100)
    print(VlnPlot(sobj, group.by = "Sample", features = feat, pt.size = 0.1) + NoLegend())
    dev.off()
}
cat(paste(dim(sobj), collapse = " x "), file=file.path(bqcdir, "dim.txt"))

print("- Applying cell QC filters ...")
print(paste("  Dim before QC:", paste(dim(sobj), collapse="x")))
if (!is.null(envs$cell_qc) && nchar(envs$cell_qc) > 0) {
    sobj = sobj %>% filter(!!rlang::parse_expr(envs$cell_qc))
}
print(paste("  Dim after QC:", paste(dim(sobj), collapse="x")))

print("- Applying gene QC filters ...")
print(paste("  Dim before QC:", paste(dim(sobj), collapse="x")))
if (is.list(envs$gene_qc)) {
    if ("min_cells" %in% names(envs$gene_qc)) {
        genes = rownames(sobj)[Matrix::rowSums(sobj) >= envs$gene_qc$min_cells]
        sobj = subset(sobj, features = genes)
    }
}
print(paste("  Dim after QC:", paste(dim(sobj), collapse="x")))

print("- Plotting QC metrics after filtering ...")
aqcdir = file.path(joboutdir, "after-qc")
dir.create(aqcdir, recursive = TRUE, showWarnings = FALSE)
for (feat in feats) {
    png(file.path(aqcdir, paste0(feat, ".png")), width = 800 + length(samples) * 15, height = 600, res = 100)
    print(VlnPlot(sobj, group.by = "Sample", features = feat, pt.size = 0.1) + NoLegend())
    dev.off()
}
cat(paste(dim(sobj), collapse = " x "), file=file.path(aqcdir, "dim.txt"))

print("- Saving results ...")
saveRDS(sobj, rdsfile)
