library(scSorter)

scsorter_db = {{envs.scsorter_db | r}}
scsorter_args = {{envs.scsorter_args | r}}
newcol = {{envs.newcol | r}}

if (is.null(scsorter_db)) { stop("`envs.scsorter_db` is not set") }

log <- get_logger()

if (startsWith(scsorter_db, "file://")) {
    scsorter_db <- sub("^file://", "", scsorter_db)
}
if (grepl("#", scsorter_db)) {
    file_path_parts <- strsplit(scsorter_db, "#")[[1]]
    scsorter_db <- file_path_parts[1]
    if (length(file_path_parts) > 1) {
        scsorter_cols <- trimws(strsplit(file_path_parts[2], ",")[[1]])
    } else {
        scsorter_cols <- NULL
    }
} else {
    scsorter_cols <- NULL
}

if (!file.exists(scsorter_db)) {
    stop(paste0("scSorter database file does not exist: ", scsorter_db))
}

log$info("Loading scSorter database ...")
if (endsWith(scsorter_db, ".rds") ||
    endsWith(scsorter_db, ".Rds") ||
    endsWith(scsorter_db, ".RDS") ||
    endsWith(scsorter_db, ".qs") ||
    endsWith(scsorter_db, ".qs2")) {
    anno <- read_obj(scsorter_db)
} else {
    anno <- read.table(scsorter_db, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

if (is.null(scsorter_cols)) {
    if (ncol(anno) < 2) {
        stop(paste0("scSorter database file must have at least 2 columns: ", scsorter_db))
    }
    if (ncol(anno) == 2) {
        colnames(anno) <- c("Type", "Marker")
    } else if (ncol(anno) >= 3) {
        colnames(anno)[1:3] <- c("Type", "Marker", "Weight")
    }
} else {
    if (length(scsorter_cols) < 2) {
        stop(paste0("scSorter database file must have at least 2 columns: ", scsorter_db))
    }
    anno <- anno[, scsorter_cols, drop = FALSE]
    if (length(scsorter_cols) == 2) {
        colnames(anno) <- c("Type", "Marker")
    } else {
        colnames(anno)[1:3] <- c("Type", "Marker", "Weight")
    }
}

log$info("Reading Seurat object...")
sobj = biopipen.utils::read_obj(sobjfile)

log$info("Running RunScSorter...")
scsorter_args$object <- sobj
scsorter_args$anno <- anno
sobj <- do_call(RunScSorter, scsorter_args)

if (!is.null(newcol)) {
    log$info(paste0("Saving cell types to new column: ", newcol, " and setting identity to it"))
    sobj@meta.data[[newcol]] <- sobj@meta.data$scSorter_celltype
    Idents(sobj) <- newcol
}

log$info("Saving Seurat object...")
biopipen.utils::save_obj(sobj, outfile)
