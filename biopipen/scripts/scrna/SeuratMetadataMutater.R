library(rlang)
library(tibble)
library(dplyr)
library(Seurat)
library(tidyseurat)
library(scplotter)
library(biopipen.utils)

srtobj <- {{in.srtobj | r}}
metafile <- {{in.metafile | r}}
outfile <- {{out.outfile | r}}
mutaters <- {{envs.mutaters | r}}
subset <- {{envs.subset | r}}

srt <- read_obj(srtobj)

if (!is.null(metafile)) {
    metadata <- srt@meta.data
    mdata <- read.table(metafile, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
    ov_cols <- intersect(colnames(metadata), colnames(mdata))
    if (length(ov_cols) > 0) {
        log_warn(paste0(
            "The following columns are already present in Seurat object and will be ignored: ",
            paste(ov_cols, collapse=', ')
        ))
    }
    srt@meta.data <- cbind(
        metadata,
        mdata[rownames(metadata), setdiff(colnames(mdata), ov_cols), drop=FALSE]
    )
}

srt <- MutateSeuratMeta(srt, mutaters)

if (!is.null(subset)) {
    srt <- filter(srt, !!parse_expr(subset))
}

save_obj(srt, outfile)
