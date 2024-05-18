source("{{biopipen_dir}}/utils/misc.R")
library(scCATCH)
library(Seurat)

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
sccatch_args = {{envs.sccatch_args | r}}
newcol = {{envs.newcol | r}}

if (!is.null(sccatch_args$marker)) {
    cellmatch = readRDS(sccatch_args$marker)
    sccatch_args$if_use_custom_marker = TRUE
}
sccatch_args$marker = cellmatch

if (is.integer(sccatch_args$use_method)) {
    sccatch_args$use_method = as.character(sccatch_args$use_method)
}

sobj = readRDS(sobjfile)

obj = createscCATCH(data = GetAssayData(sobj), cluster = as.character(Idents(sobj)))
sccatch_args$object = obj

obj = do_call(findmarkergene, sccatch_args)
obj = findcelltype(object = obj)

write.table(
    obj@celltype,
    file = file.path(dirname(outfile), "cluster2celltype.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE)

celltypes = as.list(obj@celltype$cell_type)
names(celltypes) = obj@celltype$cluster

if (length(celltypes) == 0) {
    warning("No cell types annotated from the database!")
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
}
saveRDS(sobj, outfile)
