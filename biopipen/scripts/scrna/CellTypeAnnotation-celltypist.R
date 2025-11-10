library(rlang)
library(hdf5r)
library(dplyr)
library(Seurat)
library(biopipen.utils)

sobjfile <- {{in.sobjfile | r}}
outfile <- {{out.outfile | r}}
newcol <- {{envs.newcol | r}}
cluster_ident <- {{envs.ident | r }}
merge_same_labels <- {{envs.merge | r}}
celltypist_args <- {{envs.celltypist_args | r}}
outtype <- {{envs.outtype | r }}
if (identical(outtype, "input")) {
    outtype <- tolower(tools::file_ext(outfile))  # rds, h5ad, qs/qs2
}

outdir <- dirname(outfile)
outprefix <- file.path(outdir, tools::file_path_sans_ext(basename(outfile)))

over_clustering <- celltypist_args$over_clustering %||% cluster_ident

require_package("celltypist", version = ">=1.7.1", python = celltypist_args$python)

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
ident <- NULL
if (!endsWith(sobjfile, ".h5ad")) {
    sobj <- read_obj(sobjfile)
    ident <- GetIdentityColumn(sobj)
    over_clustering <- over_clustering %||% ident

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
    # print("- {command}")
    log$debug("  {command}")
    rc <- system(command)
    if (rc != 0) {
        stop("Failed to run celltypist. Check the job.stderr file to see the error message.")
    }
}

if (outtype == "h5ad") {
    if (merge_same_labels) {
        log$warn("- Merging clusters with the same labels is not supported and is ignored for h5ad outfile ...")
    }
} else if (outtype == "rds" || outtype == "qs" || outtype == "qs2") {
    if (is.null(sobj)) {
        log$info("Reading H5AD from celltypist ...")
        sobj <- ConvertAnnDataToSeurat(
            infile = celltypist_outfile,
            outfile = NULL,
            assay = celltypist_args$assay %||% "RNA",
            ident = ident,
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

    if (celltypist_args$majority_voting) {
        prediction <- "majority_voting"

        if (!is.null(newcol)) {
            sobj@meta.data[[newcol]] <- sobj@meta.data[[prediction]]
        } else if (!isFALSE(over_clustering) && !is.null(over_clustering)) {
            # save the original over_clustering column as seurat_clusters_id
            sobj@meta.data$seurat_clusters_id <- sobj@meta.data[[over_clustering]]

            # make a map of original cluster id to new cluster id
            cluster_map <- data.frame(
                seurat_clusters_id = sobj@meta.data$seurat_clusters_id,
                seurat_clusters = sobj@meta.data[[prediction]]
                ) %>%
                group_by(seurat_clusters_id) %>%
                summarise(seurat_clusters = first(seurat_clusters), .groups = "drop") %>%
                mutate(seurat_clusters = make.unique(seurat_clusters))
            cluster_map <- split(cluster_map$seurat_clusters, cluster_map$seurat_clusters_id)
            sobj <- rename_idents(sobj, over_clustering, cluster_map)
        }
    } else if (!is.null(newcol)) {
        sobj@meta.data[[newcol]] <- sobj@meta.data[["predicted_labels"]]
    }

    if (merge_same_labels) {
        log$info("Merging clusters with the same labels ...")
        sobj <- merge_clusters_with_same_labels(sobj, newcol)
    }

    if (!is.null(ident)) {
        # restore the original identity
        Idents(sobj) <- ident
    }

    log$info("Saving the object ...")
    save_obj(sobj, outfile)
} else {
    stop(paste0("Unknown output type: ", outtype))
}
