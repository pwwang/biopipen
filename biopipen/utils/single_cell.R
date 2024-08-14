suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
try(suppressPackageStartupMessages(library(immunarch)))

#' Expand a Immunarch object into cell-level
#'
#' Here is how the data is expanded:
#' 1. Expand `$data` by Barcode (other columns are copied)
#' 2. Add sample to `Sample` column
#' 3. Add columns from `$meta`
#'
#' @param immdata Immunarch object
#' @return A data frame
#'
#' @example
#' immunarch::immdata$data$MS1 |> head()
#' #   Clones Proportion CDR3.nt            CDR3.aa V.name D.name J.name V.end D.start D.end J.start VJ.ins VD.ins DJ.ins Sequence
#' #    <dbl>      <dbl> <chr>              <chr>   <chr>  <chr>  <chr>  <int>   <int> <int>   <int>  <dbl>  <dbl>  <dbl> <lgl>
#' # 1    539    0.0634  TGTGCCAGCAGCTTACA… CASSLQ… TRBV7… TRBD2  TRBJ2…    14      18    26      29     -1      3      2 NA
#' # 2    320    0.0376  TGTGCCAGCAGCGTGTA… CASSVY… TRBV9  TRBD1  TRBJ2…    13      20    22      29     -1      6      6 NA
#' immunarch::immdata$meta |> head()
#' #   Sample  ID    Sex     Age Status Lane
#' #   <chr>   <chr> <chr> <dbl> <chr>  <chr>
#' # 1 A2-i129 C1    M        11 C      A
#' # 2 A2-i131 C2    M         9 C      A
#' # 3 A2-i133 C4    M        16 C      A
#' # 4 A2-i132 C3    F         6 C      A
#' # 5 A4-i191 C8    F        22 C      B
#' # 6 A4-i192 C9    F        24 C      B
#'
#' @export
expand_immdata <- function(immdata, cell_id = "Barcode") {
    if (!cell_id %in% colnames(immdata$data[[1]])) {
        stop(paste0("cell_id '", cell_id, "' not found in data"))
    }
    do.call(rbind, lapply(names(immdata$data), function(name) {
        # Split barcodes
        dat <- immdata$data[[name]] %>% separate_rows(!!sym(cell_id), sep = ";")
        dat$Sample <- name
        dat <- dat %>% left_join(immdata$meta, by = "Sample", suffix = c("_data", ""))

        dat
    }))
}

#' Filter expanded immdata
#'
#' @param exdata Expanded immdata
#' @param filters Filters
#' @return Filtered data
#'
#' @export
filter_expanded_immdata <- function(exdata, filters, update_clones = FALSE) {
    if (length(filters) == 0) {
        return(exdata)
    }
    out <- exdata %>% dplyr::filter(!!parse_expr(filters))
    if (update_clones) {
        out <- out %>%
            group_by(Sample, CDR3.aa) %>%
            mutate(Clones = n()) %>%
            ungroup() %>%
            group_by(Sample) %>%
            mutate(Proportion = Clones / n()) %>%
            ungroup() %>%
            arrange(Sample, desc(Clones))
    }
    out
}

#' Convert expanded immdata to Immunarch object
#'
#' @param exdata Expanded immdata
#' @param metacols Columns to be added to `$meta`
#' @return Immunarch object
#'
#' @export
immdata_from_expanded <- function(
    exdata,
    metacols = NULL,
    cell_id = "Barcode",
    update_clones = TRUE
) {
    if (is.null(metacols)) {
        metacols = setdiff(colnames(exdata), c(
            "Clones", "Proportion", "CDR3.nt", "CDR3.aa", "V.name", "D.name", "J.name",
            "V.end", "D.start", "D.end", "J.start", "VJ.ins", "VD.ins", "DJ.ins",
            "Sequence", "chain", "Barcode", "raw_clonotype_id", "ContigID", "C.name",
            "CDR1.nt", "CDR2.nt", "CDR1.aa", "CDR2.aa", "FR1.nt", "FR2.nt", "FR3.nt",
            "FR4.nt", "FR1.aa", "FR2.aa", "FR3.aa", "FR4.aa"
        ))
    }
    out <- list(meta = exdata[, metacols, drop = FALSE])
    out$meta <- out$meta[!duplicated(out$meta$Sample), , drop = FALSE]
    out$data <- lapply(
        split(
            exdata[, setdiff(colnames(exdata), metacols), drop = FALSE],
            exdata$Sample
        ),
        function(dat) {
            ncells <- nrow(dat)
            dat_cols <- setdiff(colnames(dat), c("Clones", "Proportion", cell_id))
            dat %>% group_by(CDR3.aa) %>%
                summarise(
                    Clones = ifelse(update_clones, n(), first(Clones)),
                    Proportion = ifelse(update_clones, n() / ncells, first(Proportion)),
                    !!sym(cell_id) := paste0(!!sym(cell_id), collapse = ";"),
                    !!!parse_exprs(sapply(dat_cols, function(x) paste0('first(`', x, '`)'))),
                    .groups = "drop"
                ) %>%
                arrange(desc(Clones))
        }
    )
    out
}

#' Convert Seurat object to Anndata
#'
#' @param sobjfile Seurat object file
#' @param outfile Output file
#' @param assay Assay to be used
#'
#' @export
seurat_to_anndata <- function(sobjfile, outfile, assay = NULL, log_info, tmpdir = NULL, log_indent = "") {
    library(Seurat)
    library(SeuratDisk)
    library(hdf5r)
    if (endsWith(sobjfile, ".rds") || endsWith(sobjfile, ".RDS")) {
        library(digest)

        dig <- digest::digest(sobjfile, algo = "md5")
        dig <- substr(dig, 1, 8)
        assay_name <- ifelse(is.null(assay), "", paste0("_", assay))
        tmpdir <- tmpdir %||% dirname(outfile)
        dir.create(tmpdir, showWarnings = FALSE)
        h5seurat_file <- file.path(
            tmpdir,
            paste0(
                tools::file_path_sans_ext(basename(outfile)),
                assay_name, ".", dig, ".h5seurat"
            )
        )
        if (file.exists(h5seurat_file) &&
            (file.mtime(h5seurat_file) < file.mtime(sobjfile))) {
            file.remove(h5seurat_file)
        }
        if (!file.exists(h5seurat_file)) {
            log_info("{log_indent}Reading RDS file ...")
            sobj <- readRDS(sobjfile)
            assay <- assay %||% DefaultAssay(sobj)
            # In order to convert to h5ad
            # https://github.com/satijalab/seurat/issues/8220#issuecomment-1871874649
            sobj$RNAv3 <- as(object = sobj[[assay]], Class = "Assay")
            DefaultAssay(sobj) <- "RNAv3"
            sobj$RNA <- NULL
            sobj <- RenameAssays(sobj, RNAv3 = "RNA")

            log_info("{log_indent}Saving to H5Seurat file ...")
            SaveH5Seurat(sobj, h5seurat_file)
            rm(sobj)
            gc()
            sobjfile <- h5seurat_file
        } else {
            log_info("{log_indent}Using existing H5Seurat file ...")
        }
    }

    if (!endsWith(sobjfile, ".h5seurat")) {
        stop(paste0("Unknown input file format: ",
            tools::file_ext(sobjfile),
            ". Supported formats: .rds, .RDS, .h5seurat"))
    }

    log_info("{log_indent}Converting to Anndata ...")
    Convert(sobjfile, dest = outfile, assay = assay %||% "RNA", overwrite = TRUE)

    log_info("{log_indent}Fixing categorical data ...")
    # See: https://github.com/mojaveazure/seurat-disk/issues/183
    H5.create_reference <- function(self, ...) {
        space <- self$get_space()
        do.call("[", c(list(space), list(...)))
        ref_type <- hdf5r::h5const$H5R_OBJECT
        ref_obj <- hdf5r::H5R_OBJECT$new(1, self)
        res <- .Call("R_H5Rcreate", ref_obj$ref, self$id, ".", ref_type,
                    space$id, FALSE, PACKAGE = "hdf5r")
        if (res$return_val < 0) {
            stop("Error creating object reference")
        }
        ref_obj$ref <- res$ref
        return(ref_obj)
    }

    h5ad <- H5File$new(outfile, "r+")
    cats <- names(h5ad[["obs/__categories"]])
    for (cat in cats) {
        catname <- paste0("obs/__categories/", cat)
        obsname <- paste0("obs/", cat)
        ref <- H5.create_reference(h5ad[[catname]])
        h5ad[[obsname]]$create_attr(
            attr_name = "categories",
            robj = ref,
            space = H5S$new(type = "scalar")
        )
    }
    h5ad$close()
}
