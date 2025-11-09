library(scCATCH)
library(Seurat)
library(biopipen.utils)

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
sccatch_args = {{envs.sccatch_args | r}}
newcol = {{envs.newcol | r}}
ident = {{envs.ident | r }}
merge_same_labels = {{envs.merge | r}}

log <- get_logger()

if (!is.null(sccatch_args$marker)) {
    cellmatch = read_obj(sccatch_args$marker)
    sccatch_args$if_use_custom_marker = TRUE
}
sccatch_args$marker = cellmatch

if (is.integer(sccatch_args$use_method)) {
    sccatch_args$use_method = as.character(sccatch_args$use_method)
}

log$info("Reading Seurat object...")
sobj = read_obj(sobjfile)
ident <- ident %||% GetIdentityColumn(sobj)
Idents(sobj) <- ident

log$info("Running createscCATCH ...")
obj = createscCATCH(data = GetAssayData(sobj), cluster = as.character(Idents(sobj)))
sccatch_args$object = obj

log$info("Running findmarkergene ...")
obj = do_call(findmarkergene, sccatch_args)

log$info("Running findcelltype ...")
obj = findcelltype(object = obj)

log$info("Saving the mappings ...")
write.table(
    obj@celltype,
    file = file.path(dirname(outfile), "cluster2celltype.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE)

celltypes = as.list(obj@celltype$cell_type)
names(celltypes) = obj@celltype$cluster

if (length(celltypes) == 0) {
    log$warn("- No cell types annotated from the database!")
} else {
    if (is.null(newcol)) {
        sobj <- rename_idents(sobj, ident, celltypes)
    } else {
        sobj@meta.data[[newcol]] = celltypes[as.character(Idents(sobj))]
    }

    if (merge_same_labels) {
        log$info("Merging clusters with the same labels ...")
        sobj = merge_clusters_with_same_labels(sobj, newcol)
    }
}

log$info("Saving Seurat object ...")
save_obj(sobj, outfile)
