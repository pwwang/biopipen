source("{{biopipen_dir}}/utils/misc.R")

library(rlang)
library(hdf5r)
library(dplyr)
library(Seurat)

sobjfile <- {{in.sobjfile | r}}
outfile <- {{out.outfile | r}}
newcol <- {{envs.newcol | r}}
celltypist_args <- {{envs.celltypist_args | r}}

outdir <- dirname(outfile)
outprefix <- file.path(outdir, tools::file_path_sans_ext(basename(outfile)))

if (is.null(celltypist_args$model)) {
    stop("Please specify a model for celltypist (envs.celltypist_args.model)")
} else if (!file.exists(celltypist_args$model)) {
    stop(paste0("Model file not found (envs.celltypist_args.model)"))
}
dir.create(file.path(outdir, "data", "models"), recursive = TRUE, showWarnings = FALSE)
modelfile <- file.path(outdir, "data", "models", basename(celltypist_args$model))
if (!file.exists(modelfile)) {
    file.symlink(celltypist_args$model, modelfile)
} else {
    real_modelfile <- normalizePath(Sys.readlink(modelfile))
    if (real_modelfile != normalizePath(celltypist_args$model)) {
        file.remove(modelfile)
        file.symlink(celltypist_args$model, modelfile)
    }
}

sobj <- NULL
outtype <- tolower(tools::file_ext(outfile))  # .rds, .h5ad, .h5seurat
if (!endsWith(sobjfile, ".h5ad")) {
    library(SeuratDisk)

    assay <- celltypist_args$assay
    if (endsWith(sobjfile, ".rds") || endsWith(sobjfile, ".RDS")) {
        h5s_file <- paste0(outprefix, ".h5seurat")
        if (file.exists(h5s_file) && (file.mtime(h5s_file) < file.mtime(sobjfile))) {
            file.remove(h5s_file)
        }
        if (!file.exists(h5s_file)) {
            log_info("Reading RDS file ...")
            sobj <- readRDS(sobjfile)
            assay <- assay %||% DefaultAssay(sobj)
            # In order to convert to h5ad
            # https://github.com/satijalab/seurat/issues/8220#issuecomment-1871874649
            sobj_v3 <- sobj
            sobj_v3$RNAv3 <- as(object = sobj[[assay]], Class = "Assay")
            DefaultAssay(sobj_v3) <- "RNAv3"
            sobj_v3$RNA <- NULL
            sobj_v3 <- RenameAssays(sobj_v3, RNAv3 = "RNA")

            log_info("Saving to H5Seurat file ...")
            SaveH5Seurat(sobj_v3, h5s_file)
            rm(sobj_v3)
        } else if (outtype == "rds") {
            log_info("Reading RDS file ...")
            sobj <- readRDS(sobjfile)
            assay <- assay %||% DefaultAssay(sobj)
            log_info("Using existing H5Seurat file ...")
        } else {
            log_info("Using existing H5Seurat file ...")
        }
        sobjfile <- h5s_file
    }
    if (!endsWith(sobjfile, ".h5seurat")) {
        stop(paste0("Unknown input file format: ",
            tools::file_ext(sobjfile),
            ". Supported formats: .rds, .RDS, .h5ad, .h5seurat"))
    }
    if (!endsWith(sobjfile, ".h5ad")) {  # .h5seurat
        destfile <- paste0(outprefix, ".h5ad")
        if (file.exists(destfile) && (file.mtime(destfile) < file.mtime(sobjfile))) {
            file.remove(destfile)
        }
        if (file.exists(destfile)) {
            log_info("Using existing H5AD file ...")
        } else {
            log_info("Converting to H5AD file ...")
            Convert(sobjfile, dest = destfile, assay = assay %||% "RNA")
        }
        sobjfile <- destfile
    }
}

# sobjfile h5ad ensured
# use celltypist to annotate
log_info("Annotating cell types using celltypist ...")
celltypist_script <- file.path(
    "{{biopipen_dir}}", "scripts", "scrna", "celltypist-wrapper.py"
)

if (outtype == "h5ad") {
    celltypist_outfile <- outfile
} else if (outtype == "h5seurat") {
    celltypist_outfile <- paste0(outprefix, ".celltypist.h5ad")
} else if (outtype == "rds") {
    ext <- if (is.null(sobj)) ".h5ad" else ".txt"
    celltypist_outfile <- paste0(outprefix, ".celltypist", ext)
} else {
    stop(paste0("Unknown output type: ", outtype))
}

if (file.exists(celltypist_outfile) &&
    (file.mtime(celltypist_outfile) > file.mtime(sobjfile))) {
    log_info("Using existing celltypist results ...")
} else {
    command <- paste(
        paste0("CELLTYPIST_FOLDER='", outdir, "'"),
        celltypist_args$python,
        celltypist_script,
        "-i", sobjfile,
        "-m", celltypist_args$model,
        "-o", celltypist_outfile
    )
    if (!isFALSE(celltypist_args$over_clustering) &&
        !is.null(celltypist_args$over_clustering)) {
        command <- paste(command, "-c", celltypist_args$over_clustering)
    }
    if (isTRUE(celltypist_args$majority_voting)) {
        command <- paste(command, "-v")
    }
    print("Running celltypist:")
    print(command)
    log_debug("- {command}")
    rc <- system(command)
    if (rc != 0) {
        stop("Failed to run celltypist")
    }
}

if (outtype == "h5ad") {
    # log_info("Using H5AD from celltypist as output directly ...")
    # file.rename(paste0(out_prefix, ".h5ad"), outfile)
} else if (outtype == "h5seurat") {
    log_info("Converting H5AD from celltypist to H5Seurat ...")
    # outfile is cleaned by the pipeline anyway
    Convert(
        celltypist_outfile, assay = assay %||% 'RNA', dest = outfile, overwrite = TRUE)
} else if (outtype == "rds") {
    if (is.null(sobj)) {
        log_info("Converting H5AD from celltypist to RDS ...")
        h5seurat_file <- paste0(outprefix, ".celltypist.h5seurat")
        if (file.exists(h5seurat_file) &&
            (file.mtime(h5seurat_file) > file.mtime(celltypist_outfile))) {
            log_info("- Using existing H5Seurat file ...")
        } else {
            log_info("- Converting to h5seurat ...")
            Convert(
                celltypist_outfile,
                assay = assay %||% 'RNA', dest = h5seurat_file, overwrite = TRUE)
        }
        log_info("- Converting to RDS ...")
        # Fix Missing required datasets 'levels' and 'values'
        # https://github.com/mojaveazure/seurat-disk/issues/109#issuecomment-1722394184
        f <- H5File$new(h5seurat_file, "r+")
        groups <- f$ls(recursive = TRUE)

        for (name in groups$name[grepl("categories", groups$name)]) {
            names <- strsplit(name, "/")[[1]]
            names <- c(names[1:length(names) - 1], "levels")
            new_name <- paste(names, collapse = "/")
            f[[new_name]] <- f[[name]]
        }

        for (name in groups$name[grepl("codes", groups$name)]) {
            names <- strsplit(name, "/")[[1]]
            names <- c(names[1:length(names) - 1], "values")
            new_name <- paste(names, collapse = "/")
            f[[new_name]] <- f[[name]]
            grp <- f[[new_name]]
            grp$write(args = list(1:grp$dims), value = grp$read() + 1)
        }
        f$close_all()
        # end

        sobj <- LoadH5Seurat(h5seurat_file)
        saveRDS(sobj, outfile)
    } else {
        log_info("Attaching celltypist results to Seurat object ...")

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

        if (celltypist_args$majority_voting) {
            prediction <- "majority_voting"

            if (!is.null(newcol)) {
                sobj@meta.data[[newcol]] <- sobj@meta.data[[prediction]]
            } else {
                over_clustering <- celltypist_args$over_clustering
                if (over_clustering %in% colnames(sobj@meta.data)) {
                    sobj@meta.data$seurat_clusters_id <- sobj@meta.data[[over_clustering]]
                } else {
                    over_clustering <- "over_clustering"
                }

                # make a map of original cluster id to new cluster id
                cluster_map <- data.frame(
                    seurat_clusters_id = sobj@meta.data[[over_clustering]],
                    seurat_clusters = sobj@meta.data[[prediction]]
                    ) %>%
                    group_by(seurat_clusters_id) %>%
                    summarise(seurat_clusters = first(seurat_clusters), .groups = "drop") %>%
                    mutate(seurat_clusters = make.unique(seurat_clusters))
                cluster_map <- split(cluster_map$seurat_clusters, cluster_map$seurat_clusters_id)
                if (over_clustering != "seurat_clusters") {
                    sobj@meta.data$seurat_clusters <- sobj@meta.data[[over_clustering]]
                }
                Idents(sobj) <- "seurat_clusters"
                cluster_map$object <- sobj
                log_info("Renaming clusters ...")
                sobj <- do_call(RenameIdents, cluster_map)
                sobj@meta.data$seurat_clusters <- Idents(sobj)
            }
        } else if (!is.null(newcol)) {
            sobj@meta.data[[newcol]] <- sobj@meta.data[["predicted_labels"]]
        }
        log_info("Saving Seurat object in RDS ...")
        saveRDS(sobj, outfile)
    }
} else {
    stop(paste0("Unknown output type: ", outtype))
}
