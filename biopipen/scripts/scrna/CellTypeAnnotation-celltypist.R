# CellTypeAnnotation-celltypist.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_celltypist <- function(
    sobj, ident, celltypist_args, outdir,
    h5ad_path = NULL, case_id = NULL
) {
    library(hdf5r)

    log <- get_logger()

    if (is.null(celltypist_args$model)) {
        stop("Please specify a model for celltypist (celltypist_args.model)")
    } else if (!file.exists(celltypist_args$model)) {
        stop(paste0("Model file not found (celltypist_args.model)"))
    }

    require_package(
        "celltypist2",
        version = ">=1.7.1",
        python = celltypist_args$python
    )

    over_clustering <- celltypist_args$over_clustering %||% ident

    # Use pre-converted h5ad or convert now
    if (is.null(h5ad_path)) {
        case_suffix <- if (!is.null(case_id)) paste0(".", case_id) else ""
        h5ad_path <- file.path(outdir, paste0("celltypist", case_suffix, ".h5ad"))
        log$info("Converting Seurat object to h5ad ...")
        ConvertSeuratToAnnData(
            sobj,
            outfile = h5ad_path,
            assay = celltypist_args$assay,
            log = log
        )
    }

    # Symlink the model file
    model_dir <- file.path(outdir, "data", "models")
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
    modelfile <- file.path(model_dir, basename(celltypist_args$model))
    suppressWarnings(file.remove(modelfile))
    file.symlink(normalizePath(celltypist_args$model), modelfile)

    # Determine output file
    case_suffix <- if (!is.null(case_id)) paste0("_", case_id) else ""
    celltypist_outfile <- file.path(
        outdir,
        paste0("celltypist", case_suffix, ".txt")
    )

    # Run celltypist
    biopipen_dir <- get_biopipen_dir(celltypist_args$python)
    celltypist_script <- file.path(
        biopipen_dir, "scripts", "scrna", "celltypist-wrapper.py"
    )

    if (file.exists(celltypist_outfile) &&
        (file.mtime(celltypist_outfile) > file.mtime(h5ad_path))) {
        log$warn(
            "Using existing celltypist results: {celltypist_outfile} ..."
        )
    } else {
        command <- paste(
            paste0("CELLTYPIST_FOLDER='", outdir, "'"),
            celltypist_args$python,
            celltypist_script,
            "-i", h5ad_path,
            "-m", celltypist_args$model,
            "-o", celltypist_outfile
        )
        if (!isFALSE(over_clustering) && !is.null(over_clustering)) {
            command <- paste(command, "-c", over_clustering)
        }
        if (isTRUE(celltypist_args$majority_voting)) {
            command <- paste(command, "-v")
        }
        log$info("Running celltypist:")
        print(paste0("- ", command))
        log$debug("  {command}")
        rc <- system(command)
        if (rc != 0) {
            stop(paste(
                "Failed to run celltypist.",
                "Check the job.stderr file to see the error message."
            ))
        }
    }

    # Read results
    log$info("Reading celltypist results ...")
    celltypist_out <- read.table(
        celltypist_outfile, sep = "\t", header = TRUE, row.names = 1
    )

    output_col <- ifelse(
        isTRUE(celltypist_args$majority_voting),
        "majority_voting",
        "predicted_labels"
    )

    # Extract mapping
    if (is.null(over_clustering) || isFALSE(over_clustering)) {
        # Cell-level prediction: return data frame
        annotations <- celltypist_out[, output_col, drop = FALSE]
        colnames(annotations) <- output_col
        list(cell_annotations = annotations, annotation_col = output_col)
    } else {
        # Cluster-level: build mapping from over_clustering to output_col
        combined <- data.frame(
            cluster = sobj@meta.data[
                rownames(celltypist_out), over_clustering
            ],
            celltype = celltypist_out[[output_col]],
            stringsAsFactors = FALSE
        )
        mapping <- combined %>%
            distinct(cluster, celltype)
        mapping <- stats::setNames(
            as.list(mapping$celltype),
            mapping$cluster
        )
        list(mapping = mapping)
    }
}
