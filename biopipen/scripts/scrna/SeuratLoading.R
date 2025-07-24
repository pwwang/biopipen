library(Seurat)

metafile = {{in.metafile | r}}
rdsfile = {{out.rdsfile | r}}

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

out = list()
for (i in seq_len(nrow(metadata))) {
    sample = as.character(metadata[i, "Sample", drop=T])
    path = as.character(metadata[i, "RNAData", drop=T])
    if (is.na(path) || !is.character(path) || nchar(path) == 0) {
        warning(paste0("No path found for sample: ", sample))
        next
    }
    exprs = Read10X(data.dir = path)
    out[[sample]] = CreateSeuratObject(counts=exprs)
    print(paste("Sample loaded:", sample))
}

saveRDS(out, rdsfile)
