library(hdf5r)

celltypist_args <- {{envs.celltypist_args | r}}
if (identical(outtype, "input")) {
    outtype <- tolower(tools::file_ext(outfile))  # rds, h5ad, qs/qs2
}

outdir <- dirname(outfile)
outprefix <- file.path(outdir, tools::file_path_sans_ext(basename(outfile)))

over_clustering <- celltypist_args$over_clustering %||% ident

require_package("celltypist2", version = ">=1.7.1", python = celltypist_args$python)

log <- get_logger()

if (is.null(celltypist_args$model)) {
    stop("Please specify a model for celltypist (envs.celltypist_args.model)")
} else if (!file.exists(celltypist_args$model)) {
    stop(paste0("Model file not found (envs.celltypist_args.model)"))
}
dir.create(file.path(outdir, "data", "models"), recursive = TRUE, showWarnings = FALSE)
modelfile <- file.path(outdir, "data", "models", basename(celltypist_args$model))
suppressWarnings(file.remove(modelfile))
file.symlink(normalizePath(celltypist_args$model), modelfile)

sobj <- NULL
orig_ident <- NULL
if (!endsWith(sobjfile, ".h5ad")) {
    sobj <- read_obj(sobjfile)
    orig_ident <- biopipen.utils::GetIdentityColumn(sobj)
    over_clustering <- over_clustering %||% orig_ident

    if (!isFALSE(over_clustering)) {
        destfile <- paste0(outprefix, ".", over_clustering, ".h5ad")
    } else {
        destfile <- paste0(outprefix, ".h5ad")
    }

    if (file.exists(destfile) && (file.mtime(destfile) < file.mtime(sobjfile))) {
        file.remove(destfile)
    }
    if (file.exists(destfile)) {
        log$warn("Using existing H5AD file: {destfile} ...")
    } else {
        log$info("Converting to H5AD file ...")
        ConvertSeuratToAnnData(
            sobj,
            outfile = destfile,
            assay = celltypist_args$assay,
            log = log
        )
    }
    sobjfile <- destfile
} else {
    f <- hdf5r::H5File$new(sobjfile, mode = "r")
    if ("active_ident" %in% hdf5r::h5attr_names(f)) {
        orig_ident <- hdf5r::h5attr(f, "active_ident")
    }
    f$close_all()

    over_clustering <- over_clustering %||% orig_ident
}

# sobjfile h5ad ensured
# use celltypist to annotate
log$info("Annotating cell types using celltypist ...")
# celltypist_script <- file.path(
#     "{ {biopipen_dir} }", "scripts", "scrna", "celltypist-wrapper.py"
# )
# In case this script is running in the cloud and <biopipen_dir> can not be found in there
# In stead, we use the python command, which is associated with the cloud environment,
# to get the biopipen directory
biopipen_dir <- get_biopipen_dir(celltypist_args$python)
celltypist_script <- file.path(
    biopipen_dir, "scripts", "scrna", "celltypist-wrapper.py"
)

if (outtype == "h5ad") {
    celltypist_outfile <- outfile
} else if (outtype == "rds" || outtype == "qs" || outtype == "qs2") {
    ext <- if (is.null(sobj)) ".h5ad" else ".txt"
    celltypist_outfile <- paste0(outprefix, ".celltypist", ext)
} else {
    stop(paste0("Unknown output type: ", outtype))
}

if (file.exists(celltypist_outfile) &&
    (file.mtime(celltypist_outfile) > file.mtime(sobjfile))) {
    log$warn("Using existing celltypist results: {celltypist_outfile} ...")
} else {
    command <- paste(
        paste0("CELLTYPIST_FOLDER='", outdir, "'"),
        celltypist_args$python,
        celltypist_script,
        "-i", sobjfile,
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
        stop("Failed to run celltypist. Check the job.stderr file to see the error message.")
    }
}

if (outtype == "h5ad") {
    # Create a Python script to modify the h5ad file
    modify_h5ad_script <- tempfile(fileext = ".py")

    python_code <- sprintf('
import anndata as ad
import re

# Read the h5ad file
adata = ad.read_h5ad("%s")

# Get parameters
over_clustering = %s
backup_col = %s
newcol = %s
merge = %s
output_col = "%s"

# Backup over_clustering column if needed
if over_clustering and over_clustering in adata.obs.columns and backup_col:
    adata.obs[backup_col] = adata.obs[over_clustering].astype(str)

# Get the predicted data
if output_col in adata.obs.columns:
    predicted_data = adata.obs[output_col].astype(str)

    # Apply merge transformation if needed
    if merge:
        predicted_data = predicted_data.str.replace(r"\\.\\d+$", "", regex=True)

    # Write to appropriate column
    if newcol:
        adata.obs[newcol] = predicted_data
    elif over_clustering:
        adata.obs[over_clustering] = predicted_data

# Save the modified h5ad file
adata.write_h5ad("%s")
print("H5AD file modified successfully")
',
        celltypist_outfile,
        if (is.null(over_clustering) || isFALSE(over_clustering)) "None" else paste0('"', over_clustering, '"'),
        if (is.null(backup_col)) "None" else paste0('"', backup_col, '"'),
        if (is.null(newcol)) "None" else paste0('"', newcol, '"'),
        if (merge) "True" else "False",
        ifelse(isTRUE(celltypist_args$majority_voting), "majority_voting", "predicted_labels"),
        celltypist_outfile
    )

    writeLines(python_code, modify_h5ad_script)

    # Run the Python script
    log$info("Modifying h5ad file using Python ...")
    rc <- system2(celltypist_args$python, args = modify_h5ad_script)

    # Clean up
    unlink(modify_h5ad_script)

    if (rc != 0) {
        stop("Failed to modify h5ad file using Python script")
    }

} else if (outtype == "rds" || outtype == "qs" || outtype == "qs2") {
    if (is.null(sobj)) {
        log$info("Reading H5AD from celltypist ...")
        sobj <- ConvertAnnDataToSeurat(
            infile = celltypist_outfile,
            outfile = NULL,
            assay = celltypist_args$assay %||% "RNA",
            ident = over_clustering,
            log = log
        )
    } else {
        log$info("Attaching celltypist results to Seurat object ...")

        celltypist_out <- read.table(
            celltypist_outfile, sep = "\t", header = TRUE, row.names = 1)

        sobj <- AddMetaData(
            sobj,
            celltypist_out[
                rownames(sobj@meta.data),
                setdiff(colnames(celltypist_out), colnames(sobj@meta.data)),
                drop = FALSE
            ]
        )
    }

    output_col <- ifelse(
        isTRUE(celltypist_args$majority_voting),
        "majority_voting",
        "predicted_labels"
    )

    if (is.null(over_clustering) || isFALSE(over_clustering)) {
        # keep as-is, whichever is the output_col
    } else {
        mapping <- sobj@meta.data[, c(over_clustering, output_col), drop = FALSE] %>%
            distinct(!!sym(over_clustering), !!sym(output_col))
        mapping <- stats::setNames(as.list(mapping[[output_col]]), mapping[[over_clustering]])
        sobj <- biopipen.utils::RenameSeuratIdents(
            sobj,
            ident = over_clustering,
            merge = merge,
            save_as = newcol,
            backup = backup_col,
            mapping
        )
    }

    if (!is.null(orig_ident)) {
        # restore the original identity
        Idents(sobj) <- orig_ident
    }

    log$info("Saving the object ...")
    save_obj(sobj, outfile)
} else {
    stop(paste0("Unknown output type: ", outtype))
}
