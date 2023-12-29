library(Matrix)

indir = {{in.indir | quote}}
outdir = {{out.outdir | quote}}
envs = {{envs | r}}

set.seed(envs$seed)
setwd(outdir)

logger <- function(...) {
  cat(paste(..., "\n"), file=stderr())
}

# Find the data files
mtx_file = Sys.glob(file.path(indir, "*matrix.mtx.gz"))
feat_file = c(
    Sys.glob(file.path(indir, "*genes.tsv.gz")),
    Sys.glob(file.path(indir, "*features.tsv.gz"))
)
barcode_file = Sys.glob(file.path(indir, "*barcodes.tsv.gz"))
if (length(mtx_file) == 0) {
    stop("No matrix file found in", indir)
}
if (length(mtx_file) > 1) {
    warning(paste("Multiple matrix files found in", indir, ", using the first one."))
}
if (length(feat_file) == 0) {
    stop("No feature file found in", indir)
}
if (length(feat_file) > 1) {
    warning(paste("Multiple feature files found in", indir, ", using the first one."))
}
if (length(barcode_file) == 0) {
    stop("No barcode file found in", indir)
}
if (length(barcode_file) > 1) {
    warning(paste("Multiple barcode files found in", indir, ", using the first one."))
}

mtx = readMM(mtx_file)
n_feats = nrow(mtx)
n_cells = ncol(mtx)
logger("- Dimension: Features:", n_feats, ", Cells:", n_cells)

if (envs$nfeats <= 1) {
    nfeats = as.integer(n_feats * envs$nfeats)
} else {
    nfeats = envs$nfeats
}
if (envs$ncells <= 1) {
    ncells = as.integer(n_cells * envs$ncells)
} else {
    ncells = envs$ncells
}

logger("- Identifying features to keep ...")
feats = read.table(feat_file, header=FALSE, row.names=NULL, check.names=FALSE)
feats_to_keep = c()
if (length(envs$feats_to_keep) > 0) {
    feats_to_keep = match(envs$feats_to_keep, feats[,2])
}

out_feats = unique(c(sample(1:n_feats, nfeats), feats_to_keep))
out_cells = sample(1:n_cells, ncells)
logger("- Resulting in", length(out_feats), "features and", ncells, "cells")

logger("- Subsetting matrix and saving it ...")
out_mtx = mtx[out_feats, out_cells, drop=FALSE]
out_mtx_file = file.path(outdir, "matrix.mtx")
writeMM(out_mtx, out_mtx_file)
system(paste("gzip", out_mtx_file))

logger("- Subsetting features and saving it ...")
out_feats = feats[out_feats, , drop=FALSE]
out_feats_file = gzfile(file.path(outdir, "features.tsv.gz"), "w")
write.table(out_feats, out_feats_file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
close(out_feats_file)

logger("- Subsetting barcodes and saving it ...")
barcodes = read.table(barcode_file, header=FALSE, row.names=NULL, check.names=FALSE)
out_barcodes = barcodes[out_cells, , drop=FALSE]
out_barcodes_file = gzfile(file.path(outdir, "barcodes.tsv.gz"), "w")
write.table(out_barcodes, out_barcodes_file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
close(out_barcodes_file)
