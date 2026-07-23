# CellTypeAnnotation-scbert.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_scbert <- function(
    sobj, ident, scbert_ref, scbert_model,
    scbert_label_dict, scbert_args,
    outdir, h5ad_path = NULL, case_id = NULL
) {
    log <- get_logger()

    # Validate inputs
    if (is.null(scbert_ref)) {
        stop("`scbert_ref` is required for scBERT annotation")
    }
    if (is.null(scbert_model)) {
        stop("`scbert_model` is required for scBERT annotation")
    }
    if (is.null(scbert_label_dict)) {
        stop("`scbert_label_dict` is required for scBERT annotation")
    }

    if (startsWith(scbert_model, "file://")) {
        scbert_model <- sub("^file://", "", scbert_model)
    }
    if (startsWith(scbert_label_dict, "file://")) {
        scbert_label_dict <- sub("^file://", "", scbert_label_dict)
    }
    if (!file.exists(scbert_model)) {
        stop(paste0("scBERT model checkpoint does not exist: ", scbert_model))
    }
    if (!file.exists(scbert_label_dict)) {
        stop(paste0(
            "scBERT label dictionary does not exist: ", scbert_label_dict
        ))
    }

    # Use pre-converted h5ad or convert now
    if (is.null(h5ad_path)) {
        case_suffix <- if (!is.null(case_id)) paste0(".", case_id) else ""
        h5ad_path <- file.path(
            outdir, paste0("scbert", case_suffix, ".h5ad")
        )
        log$info("Converting Seurat object to h5ad ...")
        ConvertSeuratToAnnData(
            sobj, outfile = h5ad_path,
            assay = scbert_args$assay, log = log
        )
    }

    # Output paths
    case_suffix <- if (!is.null(case_id)) paste0("_", case_id) else ""
    outfile <- file.path(
        outdir, paste0("scbert", case_suffix, ".txt")
    )

    # Run wrapper script
    biopipen_dir <- get_biopipen_dir(scbert_args$python)
    wrapper <- file.path(
        biopipen_dir, "scripts", "scrna", "scbert-wrapper.py"
    )
    command <- paste(
        scbert_args$python, wrapper,
        "-i", h5ad_path,
        "-m", scbert_model,
        "-l", scbert_label_dict,
        "-r", scbert_ref,
        "-o", outfile,
        "--bin-num", scbert_args$bin_num %||% 5,
        "--gene-num", scbert_args$gene_num %||% 16906,
        "--seed", scbert_args$seed %||% 2021
    )
    if (isTRUE(scbert_args$pos_embed)) {
        command <- paste(command, "--pos-embed")
    } else {
        command <- paste(command, "--no-pos-embed")
    }
    if (isTRUE(scbert_args$novel_type)) {
        command <- paste(
            command, "--novel-type",
            "--unassign-thres", scbert_args$unassign_thres %||% 0.5
        )
    }
    log$info("Running scBERT ...")
    rc <- system(command)
    if (rc != 0) {
        stop("Failed to run scBERT.")
    }

    # Read results (barcode-indexed)
    results <- read.table(
        outfile, sep = "\t", header = TRUE, row.names = 1
    )
    list(
        cell_annotations = results,
        annotation_col = "scbert_celltype"
    )
}
