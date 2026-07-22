# CellTypeAnnotation-scina.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_scina <- function(sobj, ident, scina_db, scina_args) {
    library(SCINA)

    log <- get_logger()

    if (is.null(scina_db)) { stop("`scina_db` is not set") }

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
    scina_args$exp <- exp
    scina_args$signatures <- signatures
    results <- do_call(SCINA, scina_args)

    # Aggregate per-cell results to cluster-level mapping (majority vote)
    log$info("Aggregating SCINA results by cluster...")
    cell_labels <- results$cell_labels
    clusters <- as.character(sobj@meta.data[[ident]])

    mapping <- list()
    for (cl in unique(clusters)) {
        cl_labels <- cell_labels[clusters == cl]
        cl_labels <- cl_labels[!is.na(cl_labels)]
        # Filter out "unknown" for majority vote
        known <- cl_labels[cl_labels != "unknown"]
        if (length(known) > 0) {
            mapping[[cl]] <- names(which.max(table(known)))
        } else {
            mapping[[cl]] <- "unknown"
        }
    }

    list(mapping = mapping)
}
