source("{{biopipen_dir}}/utils/misc.R")

library(rlang)
library(Seurat)
library(hdf5r)
library(SeuratDisk)

adfile <- {{in.adfile | r}}
outfile <- {{out.outfile | r}}
dotplot_check <- {{envs.dotplot_check | r}}
outdir <- dirname(outfile)
assay <- {{envs.assay | r}}
outtype <- {{envs.outtype | r}}

if (outtype == "rds") {
    h5seurat_file <- file.path(
        outdir,
        paste0(tools::file_path_sans_ext(basename(outfile)), ".h5seurat")
    )
} else if (outtype == "h5seurat") {
    h5seurat_file <- outfile
} else {
    stop("Unknown output file type: ", outtype)
}

if (file.exists(h5seurat_file) && (file.mtime(h5seurat_file) < file.mtime(adfile))) {
    file.remove(h5seurat_file)
}

if (!file.exists(h5seurat_file)) {
    log_info("Converting to H5Seurat file ...")
    Convert(adfile, dest = h5seurat_file, assay = assay, overwrite = TRUE)
} else {
    log_info("Using existing H5Seurat file ...")
}

if (outtype == "rds") {
    log_info("Converting to RDS file ...")
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
    if (!isFALSE(dotplot_check)) {
        log_info("Checking dotplot ...")
        dotfig <- file.path(outdir, "dotplot.png")
        if (isTRUE(dotplot_check)) {
            vobj <- FindVariableFeatures(
                sobj, selection.method = "vst", nfeatures = 2000)
            dotplot_check <- head(VariableFeatures(vobj), 10)
        } else if (is.character(dotplot_check)) {
            dotplot_check <- trimws(strsplit(dotplot_check, ",")[[1]])
        }
        png(dotfig, width = 800, height = 800, res = 70)
        p <- DotPlot(sobj, features = dotplot_check, assay = assay)
        print(p)
        dev.off()
    }
    saveRDS(sobj, outfile)
}
