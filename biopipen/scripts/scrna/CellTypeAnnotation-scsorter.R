# CellTypeAnnotation-scsorter.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_scsorter <- function(sobj, ident, scsorter_db, scsorter_args) {
    # Check if RunScSorter from the scSorter package is available
    # If not, this is the cran/scSorter package, which doesn't support RunScSorter
    # on Seurat objects
    # We need pwwang/scSorter, which is a fork of scSorter that supports Seurat objects
    if (!requireNamespace("scSorter", quietly = TRUE) ||
        !"RunScSorter" %in% getNamespaceExports("scSorter")) {
        stop(paste(
            "The scSorter package is not installed or does not support RunScSorter.",
            "Please install the pwwang/scSorter package from GitHub."
        ))
    }
    library(scSorter)

    log <- get_logger()

    if (is.null(scsorter_db)) { stop("`envs.scsorter.db` is not set") }

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
        anno <- read.table(
            scsorter_db, header = TRUE, sep = "\t", stringsAsFactors = FALSE
        )
    }

    if (is.null(scsorter_cols)) {
        if (ncol(anno) < 2) {
            stop(paste0(
                "scSorter database file must have at least 2 columns: ",
                scsorter_db
            ))
        }
        if (ncol(anno) == 2) {
            colnames(anno) <- c("Type", "Marker")
        } else if (ncol(anno) >= 3) {
            colnames(anno)[1:3] <- c("Type", "Marker", "Weight")
        }
    } else {
        if (length(scsorter_cols) < 2) {
            stop(paste0(
                "scSorter database file must have at least 2 columns: ",
                scsorter_db
            ))
        }
        anno <- anno[, scsorter_cols, drop = FALSE]
        if (length(scsorter_cols) == 2) {
            colnames(anno) <- c("Type", "Marker")
        } else {
            colnames(anno)[1:3] <- c("Type", "Marker", "Weight")
        }
    }

    log$info("Running RunScSorter...")
    # Set the active identity to the ident column
    Idents(sobj) <- ident
    scsorter_args$db <- NULL  # db is passed separately as scsorter_db
    scsorter_args$object <- sobj
    scsorter_args$anno <- anno
    sobj <- do_call(RunScSorter, scsorter_args)

    # RunScSorter stores per-cell predictions in the scSorter_celltype column;
    # aggregate to one type per cluster by majority vote
    log$info("Aggregating scSorter results by cluster...")
    mapping <- majority_vote(
        sobj@meta.data$scSorter_celltype,
        as.character(sobj@meta.data[[ident]])
    )

    list(mapping = mapping)
}
