# CellTypeAnnotation-schdeepinsight.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_schdeepinsight <- function(
    sobj, ident, schdeepinsight_ref, schdeepinsight_args,
    outdir, h5ad_path = NULL, case_id = NULL
) {
    log <- get_logger()

    if (is.null(schdeepinsight_ref)) {
        stop("`envs.schdeepinsight.ref` is not set")
    }
    if (startsWith(schdeepinsight_ref, "file://")) {
        schdeepinsight_ref <- sub("^file://", "", schdeepinsight_ref)
    }
    if (!file.exists(schdeepinsight_ref)) {
        stop(paste0(
            "scHDeepInsight reference file does not exist: ",
            schdeepinsight_ref
        ))
    }

    require_package(
        "SCHdeepinsight", python = schdeepinsight_args$python
    )

    # Use pre-converted h5ad or convert now
    if (is.null(h5ad_path)) {
        case_suffix <- if (!is.null(case_id)) paste0(".", case_id) else ""
        h5ad_path <- file.path(
            outdir, paste0("schdeepinsight", case_suffix, ".h5ad")
        )
        log$info("Converting Seurat object to h5ad ...")
        ConvertSeuratToAnnData(
            sobj, outfile = h5ad_path,
            assay = schdeepinsight_args$assay, log = log
        )
    }

    # Output paths
    case_suffix <- if (!is.null(case_id)) paste0("_", case_id) else ""
    workdir <- file.path(outdir, paste0("schdeepinsight", case_suffix))
    outfile <- file.path(
        outdir, paste0("schdeepinsight", case_suffix, ".txt")
    )

    # Run wrapper script
    biopipen_dir <- get_biopipen_dir(schdeepinsight_args$python)
    wrapper <- file.path(
        biopipen_dir, "scripts", "scrna", "schdeepinsight-wrapper.py"
    )
    command <- paste(
        schdeepinsight_args$python, wrapper,
        "-i", h5ad_path,
        "-r", schdeepinsight_ref,
        "-o", outfile,
        "-d", workdir,
        "-b", schdeepinsight_args$batch_size %||% 128,
        "--rhome", R.home()
    )
    log$info("Running scHDeepInsight ...")
    rc <- system(command)
    if (rc != 0) {
        stop("Failed to run scHDeepInsight.")
    }

    # Read results (barcode-indexed, all columns become meta.data)
    results <- read.table(
        outfile, sep = "\t", header = TRUE, row.names = 1
    )
    if (is.null(ident)) {
        list(
            cell_annotations = results,
            annotation_col = "predicted_detailed_type"
        )
    } else {
        # Aggregate per-cell results to cluster-level mapping (majority vote)
        log$info("Aggregating scHDeepInsight results by cluster...")
        mapping <- majority_vote(
            results[["predicted_detailed_type"]],
            as.character(sobj@meta.data[[ident]][rownames(results)])
        )
        list(
            mapping = mapping,
            cell_annotations = results,
            annotation_col = "predicted_detailed_type"
        )
    }
}
