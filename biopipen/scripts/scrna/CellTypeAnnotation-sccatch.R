library(scCATCH)

sccatch_args = {{envs.sccatch_args | r}}

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
if (is.null(ident)) {
    sobj@meta.data$Identity <- Idents(sobj)
    ident <- "Identity"
}

log$info("Running createscCATCH ...")
obj = createscCATCH(data = GetAssayData(sobj, assay = sccatch_args$assay), cluster = as.character(sobj@meta.data[[ident]]))
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
    sobj <- RenameSeuratIdents(sobj, celltypes, backup = backup_col, ident = ident, save_as = newcol, merge = merge)
}

log$info("Saving Seurat object ...")
save_obj(sobj, outfile)
