library(rlang)
library(tibble)
library(dplyr)
library(Seurat)

srtobj = {{in.srtobj | quote}}
metafile = {{in.metafile | r}}
mutaters = {{in.mutaters | default: {} | config: "toml" | r}}
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

srt@meta.data = mutate(metadata, !!!expr)

saveRDS(srt, rdsfile)
