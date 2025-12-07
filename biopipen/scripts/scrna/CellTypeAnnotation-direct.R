celltypes <- {{envs.cell_types | r}}
more_cell_types <- {{envs.more_cell_types | r}}

log <- biopipen.utils::get_logger()

if (is.null(celltypes) || length(celltypes) == 0) {
    log$warn("No cell types are given!")
    if (!is.null(more_cell_types) && length(more_cell_types) > 0) {
        log$warn("`envs.celltypes` is not given, won't process `envs.more_cell_types`!")
    }

    if (merge) {
        log$warn("Ignoring 'envs.merge' because no cell types are given!")
    }
    # create a symbolic link to the input file
    file.symlink(sobjfile, outfile)
} else {
    log$info("Loading Seurat object ...")
    sobj <- biopipen.utils::read_obj(sobjfile)
    ident <- ident %||% biopipen.utils::GetIdentityColumn(sobj)
    if (is.null(ident)) {
        sobj@meta.data$Identity <- Idents(sobj)
        ident <- "Identity"
    }
    idents <- sobj@meta.data[[ident]]
    if (is.factor(idents)) {
        idents <- levels(idents)
    } else {
        idents <- as.character(unique(idents))
    }
    process_celltypes <- function(ct, key = NULL) {
        if (length(ct) < length(idents)) {
            ct <- c(ct, idents[(length(ct) + 1):length(idents)])
        } else if (length(ct) > length(idents)) {
            ct <- ct[1:length(idents)]
            if (is.null(key)) {
                log$warn("The length of cell types is longer than the number of clusters!")
            } else {
                log$warn(paste0("The length of cell types for '", key, "' is longer than the number of clusters!"))
            }
        }
        for (i in seq_along(ct)) {
            if (is.na(ct[i])) next
            if (ct[i] == "-" || ct[i] == "") {
                ct[i] <- idents[i]
            } else if (ct[i] == "NA") {
                ct[i] <- NA
            }
        }
        ct <- as.list(ct)
        names(ct) <- idents
        return(ct)
    }

    if (!is.null(more_cell_types) && length(more_cell_types) > 0) {
        for (key in names(more_cell_types)) {
            ct <- more_cell_types[[key]]
            ct <- process_celltypes(ct, key)
            log$info(paste0("Adding additional cell type annotation: '", key, "' ..."))
            sobj <- RenameSeuratIdents(sobj, mapping = ct, ident = ident, save_as = key, merge = merge)
        }
    }

    celltypes <- process_celltypes(celltypes)
    sobj <- RenameSeuratIdents(sobj, mapping = celltypes, ident = ident, save_as = newcol, merge = merge, backup = backup_col)

    log$info("Saving Seurat object ...")
    biopipen.utils::save_obj(sobj, outfile)
}
