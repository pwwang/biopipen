# CellTypeAnnotation-scina.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_scina <- function(sobj, ident, scina_db, scina_args) {
    library(SCINA)

    log <- get_logger()

    if (is.null(scina_db)) { stop("`envs.scina.db` is not set") }

    if (startsWith(scina_db, "file://")) {
        scina_db <- sub("^file://", "", scina_db)
    }

    if (!file.exists(scina_db)) {
        stop(paste0("SCINA database file does not exist: ", scina_db))
    }

    log$info("Loading SCINA signature file ...")
    if (endsWith(tolower(scina_db), ".rds")) {
        signatures <- readRDS(scina_db)
    } else {
        signatures <- preprocess.signatures(scina_db)
    }

    if (!is.list(signatures) || is.null(names(signatures))) {
        stop("SCINA signatures must be a named list")
    }

    # Get expression matrix (log-normalized, genes x cells)
    log$info("Preparing expression matrix...")
    exp <- as.matrix(GetAssayData(sobj, layer = "data"))

    # Run SCINA
    log$info("Running SCINA...")
    scina_args$db <- NULL  # db is passed separately as scina_db
    scina_args$exp <- exp
    scina_args$signatures <- signatures
    results <- do_call(SCINA, scina_args)

    cell_labels <- results$cell_labels
    result <- data.frame(
        scina_celltype = unname(cell_labels),
        row.names = names(cell_labels)
    )

    if (is.null(ident)) {
        list(cell_annotations = result, annotation_col = "scina_celltype")
    } else {
        # Aggregate per-cell results to cluster-level mapping (majority vote)
        log$info("Aggregating SCINA results by cluster...")
        mapping <- majority_vote(
            cell_labels, as.character(sobj@meta.data[[ident]])
        )
        list(
            mapping = mapping,
            cell_annotations = result,
            annotation_col = "scina_celltype"
        )
    }
}
