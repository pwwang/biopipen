merge_clusters_with_same_labels <- function(sobj, newcol) {
    if (is.null(newcol)) {
        sobj@meta.data$seurat_clusters <- sub("\\.\\d+$", "", sobj@meta.data$seurat_clusters)
        Idents(sobj) <- "seurat_clusters"
    } else {
        sobj@meta.data[[newcol]] <- sub("\\.\\d+$", "", sobj@meta.data[[newcol]])
    }

    sobj
}
