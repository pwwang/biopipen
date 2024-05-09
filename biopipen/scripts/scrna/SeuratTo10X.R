library(DropletUtils)
library(Seurat)

srtobjfile = {{in.srtobj | r}}
outdir = {{out.outdir | r}}
version = {{envs.version | r}}
split_by = {{envs.split_by | r}}

srtobj = readRDS(srtobjfile)
if (!is.null(split_by)) {
    # check if split_by is a valid column
    if (is.null(srtobj[[split_by]])) {
        stop(paste0("Column ", split_by, " not found in Seurat object"))
    }

    # split Seurat object by split_by column
    objs <- SplitObject(srtobj, split.by = split_by)
    for (s in names(objs)) {
        counts <- GetAssayData(object = objs[[s]], layer = "counts")
        odir <- file.path(outdir, s)
        dir.create(odir, recursive = TRUE, showWarnings = FALSE)
        write10xCounts(odir, counts, version = version, overwrite = TRUE)
    }
} else {
    counts = GetAssayData(object = srtobj, layer = "counts")
    write10xCounts(outdir, counts, version = version, overwrite = TRUE)
}
