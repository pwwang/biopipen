library(Seurat)
library(dplyr)
library(hitype)

source("{{biopipen_dir}}/utils/misc.R")

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
tissue = {{envs.hitype_tissue | r}}
db = {{envs.hitype_db | r}}

if (is.null(db)) { stop("`envs.hitype_db` is not set") }

print("- Reading Seurat object...")
sobj = readRDS(sobjfile)

# prepare gene sets
print("- Preparing gene sets...")
if (startsWith(db, "hitypedb_")) {
    gs_list = gs_prepare(eval(as.symbol(db)), tissue)
} else {
    gs_list = gs_prepare(db, tissue)
}

# run RunHitype
print("- Running RunHitype...")
sobj = RunHitype(sobj, gs_list, threshold = 0.0, make_unique = TRUE)

print("- Renaming cell types...")
sobj$seurat_clusters_old = sobj$seurat_clusters
sobj$seurat_clusters = sobj$hitype

print("- Saving Seurat object...")
saveRDS(sobj, outfile)

print("- Saving the mappings ...")
celltypes = sobj@meta.data %>%
    group_by(seurat_clusters_old) %>%
    summarize(CellType = hitype[1]) %>%
    select(Cluster = seurat_clusters_old, CellType) %>%
    ungroup()

write.table(
    celltypes,
    file = file.path(dirname(outfile), "cluster2celltype.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
