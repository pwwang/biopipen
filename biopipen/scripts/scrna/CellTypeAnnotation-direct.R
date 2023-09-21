source("{{biopipen_dir}}/utils/misc.R")
library(Seurat)

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
celltypes = {{envs.cell_types | r}}
newcol = {{envs.newcol | r}}

if (is.null(celltypes) || length(celltypes) == 0) {
    warning("No cell types are given!")

    # create a symbolic link to the input file
    file.symlink(sobjfile, outfile)
} else {
    sobj = readRDS(sobjfile)
    idents = as.character(unique(Idents(sobj)))
    idents = idents[order(as.numeric(idents))]

    if (length(celltypes) < length(idents)) {
        celltypes = c(celltypes, idents[(length(celltypes) + 1):length(idents)])
    } else if (length(celltypes) > length(idents)) {
        celltypes = celltypes[1:length(idents)]
        warning(
            "The length of cell types is longer than the number of clusters!",
            immediate. = TRUE
        )
    }
    for (i in seq_along(celltypes)) {
        if (celltypes[i] == "-" || celltypes[i] == "") {
            celltypes[i] = idents[i]
        }
    }
    names(celltypes) = idents

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

    saveRDS(sobj, outfile)
}
