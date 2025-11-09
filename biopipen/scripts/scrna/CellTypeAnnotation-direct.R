library(Seurat)
library(rlang)
library(dplyr)
library(tidyseurat)

sobjfile <- {{in.sobjfile | r}}
outfile <- {{out.outfile | r}}
celltypes <- {{envs.cell_types | r}}
newcol <- {{envs.newcol | r}}
ident <- {{envs.ident | r }}
merge_same_labels <- {{envs.merge | r}}
more_cell_types <- {{envs.more_cell_types | r}}

log <- biopipen.utils::get_logger()

if (is.null(celltypes) || length(celltypes) == 0) {
    log$warn("No cell types are given!")
    if (!is.null(more_cell_types) && length(more_cell_types) > 0) {
        log$warn("`envs.celltypes` is not given, won't process `envs.more_cell_types`!")
    }

    if (merge_same_labels) {
        log$warn("Ignoring 'envs.merge' because no cell types are given!")
    }
    # create a symbolic link to the input file
    file.symlink(sobjfile, outfile)
} else {
    log$info("Loading Seurat object ...")
    sobj <- biopipen.utils::read_obj(sobjfile)
    ident <- ident %||% biopipen.utils::GetIdentityColumn(sobj)
    Idents(sobj) <- ident
    idents <- Idents(sobj)
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
            if (ct[i] == "-" || ct[i] == "") {
                ct[i] <- idents[i]
            }
        }
        names(ct) <- idents
        return(ct)
    }

    if (!is.null(more_cell_types) && length(more_cell_types) > 0) {
        for (key in names(more_cell_types)) {
            ct <- more_cell_types[[key]]
            ct <- process_celltypes(ct, key)
            log$info(paste0("Adding additional cell type annotation: '", key, "' ..."))
            sobj@meta.data[[key]] <- ct[as.character(Idents(sobj))]
        }
    }

    celltypes <- process_celltypes(celltypes)

    log$info("Renaming cell types ...")
    if (is.null(newcol)) {
        sobj <- rename_idents(sobj, ident, celltypes)
        log$info("Filtering clusters if NA ...")
        sobj <- filter(sobj, !!sym(ident) != "NA" & !is.na(!!sym(ident)))
    } else {
        sobj[[newcol]] <- celltypes[as.character(Idents(sobj))]
    }

    if (merge_same_labels) {
        log$info("Merging clusters with the same labels ...")
        sobj <- merge_clusters_with_same_labels(sobj, newcol)
    }

    log$info("Saving Seurat object ...")
    biopipen.utils::save_obj(sobj, outfile)
}
