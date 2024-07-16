library(Seurat)
library(dplyr)
library(hitype)

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
tissue = {{envs.hitype_tissue | r}}
db = {{envs.hitype_db | r}}
newcol = {{envs.newcol | r}}
merge_same_labels = {{envs.merge | r}}

if (is.null(db)) { stop("`envs.hitype_db` is not set") }

log_info("Reading Seurat object...")
sobj = readRDS(sobjfile)

# prepare gene sets
log_info("Preparing gene sets...")
if (startsWith(db, "hitypedb_") && !grepl(".", db, fixed = TRUE)) {
    gs_list = gs_prepare(eval(as.symbol(db)), tissue)
} else {
    gs_list = gs_prepare(db, tissue)
}

# run RunHitype
log_info("Running RunHitype...")
sobj = RunHitype(sobj, gs_list, threshold = 0.0, make_unique = TRUE)

log_info("Renaming cell types...")
hitype_levels = sobj@meta.data %>%
    select(seurat_clusters, hitype) %>%
    distinct(seurat_clusters, .keep_all = TRUE) %>%
    arrange(as.numeric(seurat_clusters)) %>%
    pull("hitype")

if (is.null(newcol)) {
    sobj$seurat_clusters_id = sobj$seurat_clusters
    sobj$seurat_clusters = factor(sobj$hitype, levels = hitype_levels)
    Idents(sobj) = "seurat_clusters"
} else {
    sobj[[newcol]] = factor(sobj$hitype, levels = hitype_levels)
}

if (merge_same_labels) {
    log_info("Merging clusters with the same labels...")
    sobj = merge_clusters_with_same_labels(sobj, newcol)
}

log_info("Saving Seurat object...")
saveRDS(sobj, outfile)

log_info("Saving the mappings ...")
if (is.null(newcol)) {
    celltypes = sobj@meta.data %>%
        group_by(seurat_clusters_id) %>%
        summarize(CellType = hitype[1]) %>%
        select(Cluster = seurat_clusters_id, CellType) %>%
        ungroup()
} else {
    celltypes = sobj@meta.data %>%
        group_by(seurat_clusters) %>%
        summarize(CellType = (!!sym(newcol))[1]) %>%
        select(Cluster = seurat_clusters, CellType) %>%
        ungroup()
}

write.table(
    celltypes,
    file = file.path(dirname(outfile), "cluster2celltype.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
