# CellTypeAnnotation-cellid.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_cellid <- function(sobj, ident, cellid_db, cellid_args) {
    library(CelliD)

    log <- get_logger()

    if (is.null(cellid_db)) {
        stop("`cellid_db` is required for CelliD annotation")
    }

    # Load marker gene list from file
    if (startsWith(cellid_db, "file://")) {
        cellid_db <- sub("^file://", "", cellid_db)
    }
    if (!file.exists(cellid_db)) {
        stop(paste0("CelliD marker gene file does not exist: ", cellid_db))
    }

    ext <- tolower(tools::file_ext(cellid_db))
    if (ext %in% c("csv", "tsv", "txt")) {
        log$info("Loading marker gene sets from CSV/TSV...")
        sep <- if (ext == "csv") "," else "\t"
        df <- read.table(
            cellid_db, header = TRUE, sep = sep, stringsAsFactors = FALSE
        )
        if (!"gene" %in% colnames(df) || !"cell_type" %in% colnames(df)) {
            stop("CSV/TSV must have 'gene' and 'cell_type' columns.")
        }
        pathways <- split(df$gene, df$cell_type)
    } else {
        log$info("Loading marker gene sets from R object file...")
        pathways <- read_obj(cellid_db)
        if (!is.list(pathways)) {
            stop("CelliD marker gene info must be a named list.")
        }
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

    list(
        cell_annotations = result,
        annotation_col = "cellid_celltype"
    )
}
