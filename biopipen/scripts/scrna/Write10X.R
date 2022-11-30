library(DropletUtils)
library(Seurat)

srtobjfile = {{in.srtobj | r}}
outdir = {{out.outdir | r}}
version = {{envs.version | r}}

srtobj = readRDS(srtobjfile)
counts = GetAssayData(object = srtobj, slot = "counts")

write10xCounts(outdir, counts, version = version, overwrite = TRUE)
