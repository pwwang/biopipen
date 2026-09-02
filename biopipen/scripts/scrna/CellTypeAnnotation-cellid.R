# CellTypeAnnotation-cellid.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_cellid <- function(sobj, ident, cellid_db, cellid_args) {
    library(CelliD)

    log <- get_logger()

    if (is.null(cellid_db)) {
        stop("`envs.cellid.db` is required for CelliD annotation")
    }

    # Load marker gene list from file
    marker_info <- load_marker_table(cellid_db)
    if (is_marker_canonical(marker_info)) {
        pathways <- markers_to_named_list(
            marker_info,
            tissue = cellid_args$tissue,
            cancer = cellid_args$cancer,
            species = cellid_args$species
        )
        cellid_args$tissue <- NULL
        cellid_args$cancer <- NULL
        cellid_args$species <- NULL
    } else if (is.data.frame(marker_info)) {
        stop("CSV/TSV must have 'gene' and 'cell_type' columns.")
    } else if (is.list(marker_info)) {
        pathways <- marker_info
    } else {
        stop("CelliD marker gene info must be a named list.")
    }

    nmcs <- cellid_args$nmcs %||% 50
    n_features <- cellid_args$n_features %||% 200
    dims <- cellid_args$dims %||% seq(nmcs)
    min_size <- cellid_args$min_size %||% 10
    log_trans <- cellid_args$log_trans %||% TRUE
    p_adjust <- cellid_args$p_adjust %||% TRUE

    # Run MCA
    log$info("Running MCA with {nmcs} components ...")
    sobj <- RunMCA(sobj, nmcs = nmcs)

    # Run per-cell hypergeometric test
    log$info(
        "Running per-cell hypergeometric test against {length(pathways)} gene sets ..."
    )
    enrichment <- RunCellHGT(
        X = sobj,
        pathways = pathways,
        reduction = "mca",
        n.features = n_features,
        dims = dims,
        minSize = min_size,
        log.trans = log_trans,
        p.adjust = p_adjust
    )

    # Convert enrichment matrix to cell type predictions (argmax per cell)
    # enrichment is: cells x pathways (cell types)
    enrichment <- as.matrix(enrichment)
    cell_type_names <- colnames(enrichment)
    pred_idx <- apply(enrichment, 1, which.max)
    predicted <- cell_type_names[pred_idx]
    names(predicted) <- rownames(enrichment)

    # Build cell annotations data frame
    result <- data.frame(
        cellid_celltype = predicted,
        row.names = names(predicted)
    )

    if (is.null(ident)) {
        list(cell_annotations = result, annotation_col = "cellid_celltype")
    } else {
        # Aggregate per-cell results to cluster-level mapping (majority vote)
        log$info("Aggregating CelliD results by cluster...")
        mapping <- majority_vote(
            predicted, as.character(sobj@meta.data[[ident]])
        )
        list(
            mapping = mapping,
            cell_annotations = result,
            annotation_col = "cellid_celltype"
        )
    }
}
