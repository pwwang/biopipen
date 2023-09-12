source("{{biopipen_dir}}/utils/mutate_helpers.R")
library(rlang)
library(tibble)
library(dplyr)
library(Seurat)

srtobj = {{in.srtobj | quote}}
metafile = {{in.metafile | r}}
mutaters = {{envs.mutaters | r}}
rdsfile = {{out.rdsfile | quote}}

srt = readRDS(srtobj)
metadata = srt@meta.data

if (!is.null(metafile)) {
    mdata = read.table(metafile, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
    metadata = cbind(metadata, mdata[rownames(metadata),,drop=FALSE])
}

expr = list()
for (key in names(mutaters)) {
    expr[[key]] = parse_expr(mutaters[[key]])
}

if (!is.null(expr) && length(expr) > 0) {
    srt@meta.data = mutate(metadata, !!!expr)
} else {
    srt@meta.data = metadata
}

saveRDS(srt, rdsfile)
