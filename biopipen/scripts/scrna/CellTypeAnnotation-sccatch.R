library(scCATCH)
library(Seurat)

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
sccatch_args = {{envs.sccatch_args | r}}
newcol = {{envs.newcol | r}}
merge_same_labels = {{envs.merge | r}}

if (!is.null(sccatch_args$marker)) {
    cellmatch = readRDS(sccatch_args$marker)
    sccatch_args$if_use_custom_marker = TRUE
}
sccatch_args$marker = cellmatch

if (is.integer(sccatch_args$use_method)) {
    sccatch_args$use_method = as.character(sccatch_args$use_method)
}

log_info("Reading Seurat object...")
sobj = readRDS(sobjfile)

log_info("Running createscCATCH ...")
obj = createscCATCH(data = GetAssayData(sobj), cluster = as.character(Idents(sobj)))
sccatch_args$object = obj

log_info("Running findmarkergene ...")
obj = do_call(findmarkergene, sccatch_args)

log_info("Running findcelltype ...")
obj = findcelltype(object = obj)

log_info("Saving the mappings ...")
write.table(
    obj@celltype,
    file = file.path(dirname(outfile), "cluster2celltype.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE)

celltypes = as.list(obj@celltype$cell_type)
names(celltypes) = obj@celltype$cluster

if (length(celltypes) == 0) {
    log_warn("- No cell types annotated from the database!")
} else {
    if (is.null(newcol)) {
        sobj$seurat_clusters_id = Idents(sobj)
        celltypes$object = sobj
        sobj = do_call(RenameIdents, celltypes)
        sobj$seurat_clusters = Idents(sobj)
    } else {
        celltypes$object = sobj
        sobj = do_call(RenameIdents, celltypes)
        sobj[[newcol]] = Idents(sobj)
        Idents(sobj) = "seurat_clusters"
    }

    if (merge_same_labels) {
        log_info("Merging clusters with the same labels ...")
        sobj = merge_clusters_with_same_labels(sobj, newcol)
    }
}

log_info("Saving Seurat object ...")
saveRDS(sobj, outfile)
