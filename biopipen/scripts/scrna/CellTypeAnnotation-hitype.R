library(rlang)
library(dplyr)
library(hitype)

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
tissue = {{envs.hitype_tissue | r}}
db = {{envs.hitype_db | r}}
newcol = {{envs.newcol | r}}
ident = {{envs.ident | r }}
merge_same_labels = {{envs.merge | r}}

if (is.null(db)) { stop("`envs.hitype_db` is not set") }

log <- get_logger()

log$info("Reading Seurat object...")
sobj = biopipen.utils::read_obj(sobjfile)
ident <- ident %||% biopipen.utils::GetIdentityColumn(sobj)
Idents(sobj) <- ident

# prepare gene sets
log$info("Preparing gene sets...")
if (startsWith(db, "hitypedb_") && !grepl(".", db, fixed = TRUE)) {
    gs_list = gs_prepare(eval(as.symbol(db)), tissue)
} else {
    gs_list = gs_prepare(db, tissue)
}

# run RunHitype
log$info("Running RunHitype...")
sobj = RunHitype(sobj, gs_list, threshold = 0.0, make_unique = TRUE)

log$info("Renaming cell types...")
hitype_labels <- sobj@meta.data %>%
    distinct(!!sym(ident), hitype)
hitype_labels <- split(hitype_labels$hitype, hitype_labels[[ident]])

if (is.null(newcol)) {
    sobj <- rename_idents(sobj, ident, hitype_labels)
} else {
    sobj[[newcol]] = sobj$hitype
}

if (merge_same_labels) {
    log$info("Merging clusters with the same labels...")
    sobj = merge_clusters_with_same_labels(sobj, newcol)
}

log$info("Saving Seurat object...")
biopipen.utils::save_obj(sobj, outfile)

log$info("Saving the mappings ...")
if (is.null(newcol)) {
    celltypes = sobj@meta.data %>%
        group_by(!!sym(backup_col)) %>%
        summarize(CellType = hitype[1]) %>%
        select(Cluster = !!sym(backup_col), CellType) %>%
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
