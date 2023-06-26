source("{{biopipen_dir}}/utils/misc.R")
library(Seurat)

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
celltypes = {{envs.cell_types | r}}

sobj = readRDS(sobjfile)
idents = as.character(unique(Idents(sobj)))
idents = idents[order(as.numeric(idents))]

if (!is.list(celltypes) && length(celltypes) > 0) {
    celltypes = celltypes[1:length(idents)]
    celltypes = as.list(celltypes)
    names(celltypes) = idents
}

if (length(celltypes) == 0) {
    warning("No cell types are given!")
} else {
    celltypes$object = sobj
    sobj = do_call(RenameIdents, celltypes)
    sobj$seurat_clusters = Idents(sobj)

}

saveRDS(sobj, outfile)
