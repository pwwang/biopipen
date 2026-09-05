# CellTypeAnnotation-cellassign.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_cellassign <- function(sobj, ident, cellassign_db, cellassign_args) {
    python <- cellassign_args$python %||% Sys.which("python")
    if (python == "") {
        stop("Python executable not found. Please specify `envs.cellassign.python`.")
    }
    # load the right Python environment with tensorflow installed
    Sys.setenv(RETICULATE_PYTHON = python)

    library(cellassign)

    log <- get_logger()

    if (is.null(cellassign_db)) {
        stop("`envs.cellassign.db` is required for cellassign annotation")
    }

    # Load marker gene info
    marker_gene_info <- load_marker_table(cellassign_db)
    if (!is_marker_canonical(marker_gene_info)) {
        # Native signature formats have no tissue/cancer/species columns
        stop_on_filtering_native_db(
            cellassign_args$tissue,
            cellassign_args$cancer,
            cellassign_args$species
        )
    }
    if (is_marker_canonical(marker_gene_info)) {
        marker_gene_info <- markers_to_named_list(
            marker_gene_info,
            tissue = cellassign_args$tissue,
            cancer = cellassign_args$cancer,
            species = cellassign_args$species
        )
    } else if (is.data.frame(marker_gene_info)) {
        stop("CSV/TSV must have 'gene' and 'cell_type' columns.")
    } else if (!(is.list(marker_gene_info) || is.matrix(marker_gene_info))) {
        stop("Marker gene info must be a named list or binary matrix.")
    }
    # The filter envs are consumed by the marker conversion above; never
    # forward them to cellassign() which has no such arguments
    cellassign_args$tissue <- NULL
    cellassign_args$cancer <- NULL
    cellassign_args$species <- NULL

    # Get raw counts
    assay <- cellassign_args$assay %||% DefaultAssay(sobj)
    log$info("Extracting raw counts from assay: {assay}")
    # genes x cells matrix
    counts <- as.matrix(GetAssayData(sobj, assay = assay, layer = "counts"))
    library_size <- colSums(counts)
    size_factors <- library_size / median(library_size)

    # Filter to marker genes present in data
    if (is.list(marker_gene_info)) {
        all_markers <- unique(unlist(marker_gene_info))
    } else if (is.matrix(marker_gene_info)) {
        all_markers <- rownames(marker_gene_info)
    } else {
        stop("marker_gene_info must be a named list or binary matrix.")
    }
    common_genes <- intersect(all_markers, rownames(counts))
    if (length(common_genes) == 0) {
        stop("None of the marker genes are found in the expression data.")
    }
    log$info(
        "Found {length(common_genes)} / {length(all_markers)} marker genes in data"
    )

    cellassign_args$s <- size_factors
    # Filter counts and marker_gene_info to common genes
    counts <- counts[common_genes, , drop = FALSE]
    if (is.list(marker_gene_info)) {
        marker_gene_info <- lapply(
            marker_gene_info,
            function(gs) intersect(gs, common_genes)
        )
        marker_gene_info <- marker_gene_info[lengths(marker_gene_info) > 0]
    } else {
        marker_gene_info <- marker_gene_info[common_genes, , drop = FALSE]
    }

    # Extract extra args for cellassign()
    extra_args <- cellassign_args
    extra_args$assay <- NULL  # already handled
    extra_args$db <- NULL  # db is passed separately as cellassign_db
    extra_args$python <- NULL  # already handled above

    # Run cellassign
    log$info("Running cellassign...")
    fit <- do_call(cellassign::cellassign, c(
        list(
            exprs_obj = t(counts),
            marker_gene_info = marker_gene_info
        ),
        extra_args
    ))

    # Build cell annotations
    cell_types <- fit$cell_type
    result <- data.frame(
        cellassign_celltype = cell_types,
        row.names = names(cell_types)
    )

    if (is.null(ident)) {
        list(cell_annotations = result, annotation_col = "cellassign_celltype")
    } else {
        # Aggregate per-cell results to cluster-level mapping (majority vote)
        log$info("Aggregating cellassign results by cluster...")
        mapping <- majority_vote(
            cell_types, as.character(sobj@meta.data[[ident]])
        )
        list(
            mapping = mapping,
            cell_annotations = result,
            annotation_col = "cellassign_celltype"
        )
    }
}
