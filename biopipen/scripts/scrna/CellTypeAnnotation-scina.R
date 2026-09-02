# CellTypeAnnotation-scina.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_scina <- function(sobj, ident, scina_db, scina_args) {
    library(SCINA)

    log <- get_logger()

    if (is.null(scina_db)) { stop("`envs.scina.db` is not set") }

    log$info("Loading SCINA signature file ...")
    mt <- load_marker_table(scina_db)
    if (is_marker_canonical(mt)) {
        signatures <- markers_to_named_list(
            mt,
            tissue = scina_args$tissue,
            cancer = scina_args$cancer,
            species = scina_args$species
        )
        scina_args$tissue <- NULL
        scina_args$cancer <- NULL
        scina_args$species <- NULL
    } else if (is.data.frame(mt)) {
        # Native per-cell-type-column CSV/TSV (one column per cell type)
        signatures <- lapply(mt, function(x) x[!is.na(x) & x != ""])
    } else if (is.list(mt)) {
        signatures <- mt  # RDS named list
    } else {
        stop("Cannot recognize the SCINA signature file format.")
    }

    if (!is.list(signatures) || is.null(names(signatures))) {
        stop("SCINA signatures must be a named list")
    }

    # Get expression matrix (log-normalized, genes x cells)
    log$info("Preparing expression matrix...")
    exp <- as.matrix(GetAssayData(sobj, layer = "data"))

    # Drop signature genes not in the expression matrix, and simulate SCINA's
    # rm_overlap removal, so signatures emptied by either are dropped here
    # with a clear warning (SCINA chokes on an empty signature with a cryptic
    # chol() error).
    log$info("Filtering signatures against the expression matrix...")
    signatures <- lapply(
        signatures,
        function(x) unique(x[x %in% row.names(exp)])
    )
    rm_overlap <- scina_args$rm_overlap %||% formals(SCINA)$rm_overlap %||% 1
    if (isTRUE(rm_overlap) || rm_overlap == 1) {
        counts <- table(unlist(signatures))
        signatures <- lapply(
            signatures,
            function(x) x[x %in% names(counts[counts == 1])]
        )
    }
    signatures <- lapply(
        signatures,
        function(x) x[apply(exp[x, , drop = F], 1, sd) > 0]
    )
    n_genes <- lengths(signatures)
    if (any(n_genes == 0)) {
        log$warn(sprintf(
            paste(
                "Dropping %d signature(s) with no genes in the expression matrix:",
                "%s"
            ),
            sum(n_genes == 0),
            paste(names(n_genes)[n_genes == 0], collapse = ", ")
        ))
        signatures <- signatures[n_genes > 0]
    }
    if (length(signatures) == 0) {
        stop("No SCINA signatures have genes in the expression matrix.")
    }

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
