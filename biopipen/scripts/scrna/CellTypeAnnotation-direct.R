library(Seurat)

sobjfile <- {{in.sobjfile | r}}
outfile <- {{out.outfile | r}}
celltypes <- {{envs.cell_types | r}}
newcol <- {{envs.newcol | r}}
merge_same_labels <- {{envs.merge | r}}

log <- biopipen.utils::get_logger()

if (is.null(celltypes) || length(celltypes) == 0) {
    log$warn("No cell types are given!")

    if (merge_same_labels) {
        log$warn("Ignoring 'envs.merge' because no cell types are given!")
    }
    # create a symbolic link to the input file
    file.symlink(sobjfile, outfile)
} else {
    log$info("Loading Seurat object ...")
    sobj <- biopipen.utils::read_obj(sobjfile)
    idents <- Idents(sobj)
    if (is.factor(idents)) {
        idents <- levels(idents)
    } else {
        idents <- as.character(unique(idents))
    }

    if (length(celltypes) < length(idents)) {
        celltypes <- c(celltypes, idents[(length(celltypes) + 1):length(idents)])
    } else if (length(celltypes) > length(idents)) {
        celltypes <- celltypes[1:length(idents)]
        log$warn("The length of cell types is longer than the number of clusters!")
    }
    for (i in seq_along(celltypes)) {
        if (celltypes[i] == "-" || celltypes[i] == "") {
            celltypes[i] <- idents[i]
        }
    }
    names(celltypes) <- idents

    log$info("Renaming cell types ...")
    if (is.null(newcol)) {
        has_na <- "NA" %in% unlist(celltypes) || anyNA(unlist(celltypes))
        sobj$seurat_clusters_id <- Idents(sobj)
        celltypes$object <- sobj
        sobj <- do_call(RenameIdents, celltypes)
        sobj$seurat_clusters <- Idents(sobj)
        if (has_na) {
            log$info("Filtering clusters if NA ...")
            sobj <- subset(
                sobj,
                subset = seurat_clusters != "NA" & !is.na(seurat_clusters)
            )
        }
    } else {
        celltypes$object <- sobj
        sobj <- do_call(RenameIdents, celltypes)
        sobj[[newcol]] <- Idents(sobj)
        Idents(sobj) <- "seurat_clusters"
    }

    if (merge_same_labels) {
        log$info("Merging clusters with the same labels ...")
        sobj <- merge_clusters_with_same_labels(sobj, newcol)
    }

    biopipen.utils::save_obj(sobj, outfile)
}
