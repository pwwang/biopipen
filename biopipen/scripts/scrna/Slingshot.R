library(rlang)
library(Seurat)
library(slingshot)
library(biopipen.utils)

sobjfile <- {{in.sobjfile | r}}
outfile <- {{out.outfile | r}}
group_by <- {{envs.group_by | r}}
reduction <- {{envs.reduction | r}}
dims <- {{envs.dims | r}}
start <- {{envs.start | r}}
end <- {{envs.end | r}}
prefix <- {{envs.prefix | r}}
reverse <- {{envs.reverse | r}}
align_start <- {{envs.align_start | r}}
seed <- {{envs.seed | r}}

set.seed(seed)

log <- get_logger()

log$info("Reading Seurat object ...")
srt <- read_obj(sobjfile)
group_by <- group_by %||% biopipen.utils::GetIdentityColumn(srt)

if (is.null(group_by) || !group_by %in% colnames(srt@meta.data)) {
    stop(paste("Grouping column", group_by, "not found in the Seurat object"))
}

reduction <- reduction %||% DefaultDimReduc(srt)
dims <- biopipen.utils:::.expand_number(dims)

if (is.null(prefix)) {
    prefix <- ""
} else {
    prefix <- paste0(prefix, "_")
}

log$info("Filtering cells in NA group_by ...")
srt_sub <- srt[, !is.na(srt[[group_by, drop = TRUE]])]

log$info("Running Slingshot ...")
sl <- slingshot(
    data = as.data.frame(srt_sub[[reduction]]@cell.embeddings[, dims]),
    clusterLabels = as.character(srt_sub[[group_by, drop = TRUE]]),
    start.clus = start, end.clus = end
)

df <- as.data.frame(slingPseudotime(sl))
colnames(df) <- paste0(prefix, colnames(df))
if (isTRUE(reverse)) {
    if (isTRUE(align_start)) {
        df <- apply(df, 2, function(x) max(x, na.rm = TRUE) - x)
    } else {
        df <- max(df, na.rm = TRUE) - df
    }
}

srt <- AddMetaData(srt, metadata = df)
srt <- AddMetaData(srt, metadata = slingBranchID(sl), col.name = paste0(prefix, "BranchID"))

srt <- AddSeuratCommand(srt, "Slingshot", "slingshot(...)")

log$info("Saving Seurat object ...")
save_obj(srt, outfile)
