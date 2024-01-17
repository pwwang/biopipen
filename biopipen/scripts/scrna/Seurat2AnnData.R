source("{{biopipen_dir}}/utils/misc.R")

library(rlang)
library(Seurat)
library(SeuratDisk)

sobjfile <- {{in.sobjfile | r}}
outfile <- {{out.outfile | r}}
outdir <- dirname(outfile)
assay <- {{envs.assay | r}}

if (endsWith(sobjfile, ".rds") || endsWith(sobjfile, ".RDS")) {
    assay_name <- ifelse(is.null(assay), "", paste0("_", assay))
    h5seurat_file <- file.path(
        outdir,
        paste0(tools::file_path_sans_ext(basename(outfile)), assay_name, ".h5seurat")
    )
    if (file.exists(h5seurat_file) &&
        (file.mtime(h5seurat_file) < file.mtime(sobjfile))) {
        file.remove(h5seurat_file)
    }
    if (!file.exists(h5seurat_file)) {
        log_info("Reading RDS file ...")
        sobj <- readRDS(sobjfile)
        assay <- assay %||% DefaultAssay(sobj)
        # In order to convert to h5ad
        # https://github.com/satijalab/seurat/issues/8220#issuecomment-1871874649
        sobj$RNAv3 <- as(object = sobj[[assay]], Class = "Assay")
        DefaultAssay(sobj) <- "RNAv3"
        sobj$RNA <- NULL
        sobj <- RenameAssays(sobj, RNAv3 = "RNA")

        log_info("Saving to H5Seurat file ...")
        SaveH5Seurat(sobj, h5seurat_file)
        rm(sobj)
        sobjfile <- h5seurat_file
    } else {
        log_info("Using existing H5Seurat file ...")
    }
}

if (!endsWith(sobjfile, ".h5seurat")) {
    stop(paste0("Unknown input file format: ",
        tools::file_ext(sobjfile),
        ". Supported formats: .rds, .RDS, .h5seurat"))
}

Convert(sobjfile, dest = outfile, assay = assay %||% "RNA", overwrite = TRUE)
