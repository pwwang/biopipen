library(dplyr)
library(HGNChelper)
library(Seurat)
source("{{biopipen_dir}}/scripts/scrna/sctype.R")

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
celltypes = {{envs.cell_types | r}}

sobj = readRDS(sobjfile)
idents = as.character(unique(Idents(sobj)))
idents = idents[order(as.numeric(idents))]

if (!is.list(celltypes)) {
    celltypes = celltypes[1:length(idents)]
    celltypes = as.list(celltypes)
    names(celltypes) = idents
}
celltypes$object = sobj

sobj = do.call(RenameIdents, celltypes)
sobj$seurat_clusters = Idents(sobj)

saveRDS(sobj, outfile)
