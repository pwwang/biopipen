# CellTypeAnnotation-singler.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_singler <- function(sobj, ident, singler_db, singler_args) {
    library(SingleR)

    log <- get_logger()

    if (is.null(singler_db)) { stop("`envs.singler.db` is not set") }

    if (startsWith(singler_db, "file://")) {
        singler_db <- sub("^file://", "", singler_db)
    }

    if (!file.exists(singler_db)) {
        stop(paste0("SingleR database file does not exist: ", singler_db))
    }

    # Detect which SingleR API is available:
    # - Bioconductor (LTLA): SingleR(test, ref, labels, clusters, ...)
    # - CRAN (dviraran):    SingleR(sc_data, ref_data, types, clusters, ...)
    is_bioc <- "test" %in% names(formals(SingleR))
    log$info(
        "Using SingleR {ifelse(is_bioc, 'Bioconductor', 'CRAN')} API"
    )

    # Load reference (supports RDS, qs, qs2 via biopipen.utils)
    log$info("Loading SingleR reference ...")
    ref <- read_obj(singler_db)

    # Resolve label column
    label_col <- singler_args$label
    singler_args$label <- NULL
    singler_args$db <- NULL  # db is passed separately as singler_db

    # Prepare reference data and labels based on API version
    if (is_bioc) {
        # Bioconductor API: ref is a SummarizedExperiment, labels is a vector
        if (inherits(ref, "Seurat")) {
            log$info("Converting Seurat reference to SummarizedExperiment ...")
            ref <- as.SingleCellExperiment(ref)
        }

        if (!is.null(label_col)) {
            if (is(ref, "SummarizedExperiment")) {
                labels <- SummarizedExperiment::colData(ref)[[label_col]]
            } else {
                stop(paste0("Label column '", label_col, "' not found"))
            }
        } else if (is(ref, "SummarizedExperiment")) {
            labels <- NULL
            for (col in c(
                "label.main", "label.fine", "label.ont", "label"
            )) {
                if (col %in% colnames(SummarizedExperiment::colData(ref))) {
                    labels <- SummarizedExperiment::colData(ref)[[col]]
                    log$info("Auto-detected label column: {col}")
                    break
                }
            }
        } else {
            stop(paste(
                "Reference must be a Seurat or SummarizedExperiment object."
            ))
        }

        if (is.null(labels)) {
            stop(paste(
                "Cannot determine labels from reference.",
                "Set `label` in `envs.singler.label`."
            ))
        }

        # Bioconductor SingleR call
        log$info("Preparing expression matrix ...")
        exp <- as.matrix(GetAssayData(sobj, layer = "data"))

        clusters <- as.character(sobj@meta.data[[ident]])
        log$info(
            "Running SingleR with {length(unique(clusters))} clusters ..."
        )

        singler_args$test <- exp
        singler_args$ref <- ref
        singler_args$labels <- labels
        singler_args$clusters <- clusters

        results <- do_call(SingleR, singler_args)

        # Build mapping (prefer pruned.labels, fall back to labels)
        mapping <- as.list(results$pruned.labels)
        names(mapping) <- rownames(results)
        na_mask <- is.na(mapping) | mapping == "NA"
        if (any(na_mask)) {
            mapping[na_mask] <- as.list(results$labels[na_mask])
        }

    } else {
        # CRAN API: ref_data is a matrix, types is a vector
        if (inherits(ref, "Seurat")) {
            log$info("Extracting data from Seurat reference ...")
            ref_data <- as.matrix(GetAssayData(ref, layer = "data"))
            meta <- ref@meta.data
        } else if (is(ref, "SummarizedExperiment")) {
            assay_name <- "logcounts"
            if (!assay_name %in% names(SummarizedExperiment::assays(ref))) {
                assay_name <- names(SummarizedExperiment::assays(ref))[1]
            }
            ref_data <- as.matrix(
                SummarizedExperiment::assay(ref, assay_name)
            )
            meta <- SummarizedExperiment::colData(ref)
        } else {
            stop(paste(
                "Reference must be a Seurat or SummarizedExperiment object."
            ))
        }

        if (!is.null(label_col)) {
            if (!label_col %in% colnames(meta)) {
                stop(paste0(
                    "Label column '", label_col, "' not found in reference"
                ))
            }
            types <- meta[[label_col]]
        } else {
            types <- NULL
            for (col in c(
                "label.main", "label.fine", "label.ont", "label"
            )) {
                if (col %in% colnames(meta)) {
                    types <- meta[[col]]
                    log$info("Auto-detected label column: {col}")
                    break
                }
            }
        }

        if (is.null(types)) {
            stop(paste(
                "Cannot determine cell types from reference.",
                "Set `label` in `envs.singler.label`."
            ))
        }

        # CRAN SingleR call
        log$info("Preparing expression matrix ...")
        sc_data <- as.matrix(GetAssayData(sobj, layer = "data"))

        clusters <- as.factor(sobj@meta.data[[ident]])
        log$info(
            "Running SingleR with {length(levels(clusters))} clusters ..."
        )

        singler_args$method <- "cluster"
        singler_args$sc_data <- sc_data
        singler_args$ref_data <- ref_data
        singler_args$types <- types
        singler_args$clusters <- clusters

        results <- do_call(SingleR, singler_args)

        # Build mapping from labels
        mapping <- as.list(results$labels)
    }

    list(mapping = mapping)
}
