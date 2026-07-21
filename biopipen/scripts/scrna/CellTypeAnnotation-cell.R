ident <- {{envs.ident | r}}
celltypes <- {{envs.cell_types | r}}

if (!is.character(celltypes)) {
    stop("`envs.cell_types` must be path to a TSV file with cell-level annotations.")
}

if (startsWith(celltypes, "file://")) {
    celltypes <- sub("^file://", "", celltypes)
}

log <- biopipen.utils::get_logger()

log$info(paste0("Loading cell type file: ", celltypes))

file_path_parts <- strsplit(celltypes, "#")[[1]]
celltypes_file <- file_path_parts[1]
if (length(file_path_parts) > 1) {
    celltypes_cols <- trimws(strsplit(file_path_parts[2], ",")[[1]])
} else {
    celltypes_cols <- c("1", "2")
}

if (!file.exists(celltypes_file)) {
    stop(paste0("Cell type file does not exist: ", celltypes_file))
}

celltypes <- read.table(celltypes_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
celltypes_cols <- sapply(
    celltypes_cols,
    function(col) {
        if (grepl("^[+-]?[0-9]+$", col)) {
            col <- as.integer(col)
            if (col < 0) {
                col <- ncol(celltypes) + col + 1
            }
            if (col < 1 || col > ncol(celltypes)) {
                stop(paste0("Column index '", col, "' is out of bounds for cell type file: ", celltypes_file))
            }
            col <- colnames(celltypes)[col]
        }
        if (!col %in% colnames(celltypes)) {
            stop(paste0("Column '", col, "' not found in cell type file: ", celltypes_file))
        }
        return(col)
    }
)
celltypes <- celltypes[, celltypes_cols, drop = FALSE]
cell_id_col <- celltypes_cols[1]
cell_type_cols <- celltypes_cols[-1]
ident <- ident %||% cell_type_cols[1]
if (!is.null(ident)) {
    if (!ident %in% cell_type_cols) {
        stop(paste0("Identity column '", ident, "' not found in cell type file: ", celltypes_file))
    }
}

log$info("Reading Seurat object...")
sobj <- biopipen.utils::read_obj(sobjfile)

existing_cells <- intersect(rownames(sobj@meta.data), celltypes[[cell_id_col]])
if (length(existing_cells) == 0) {
    stop(
        paste0(
            "No cells in the Seurat object match the cell IDs in the cell type file. \n",
            "Cell IDs in the Seurat object look like: ", paste(head(rownames(sobj@meta.data)), collapse = ", "), " ... \n",
            "Cell IDs in the cell type file look like: ", paste(head(celltypes[[cell_id_col]]), collapse = ", "), " ... \n"
        )
    )
}

existing_cols <- intersect(cell_type_cols, colnames(sobj@meta.data))
if (length(existing_cols) > 0) {
    log$warn(paste0("Column(s) '", paste(existing_cols, collapse = ", "), "' already exist(s) in Seurat object metadata, will be overwritten."))
    # remove existing columns from sobj@meta.data
    sobj@meta.data[, existing_cols] <- NULL
}


missing_cells <- setdiff(rownames(sobj@meta.data), celltypes[[cell_id_col]])
if (length(missing_cells) > 0) {
    log$warn(
        paste0("The following cells are missing in the cell type file: ",
        paste(head(missing_cells), collapse = ", "),
        if (length(missing_cells) > 6) paste0(", ... (", length(missing_cells) - 6, " more)"), ".")
    )
}
# Add NAs for missing cells
missing_celltypes <- data.frame(matrix(NA, nrow = length(missing_cells), ncol = ncol(celltypes)))
colnames(missing_celltypes) <- colnames(celltypes)
missing_celltypes[[cell_id_col]] <- missing_cells
celltypes <- rbind(celltypes, missing_celltypes)

# Reorder celltypes to match the order of cells in the Seurat object
celltypes <- celltypes[match(rownames(sobj@meta.data), celltypes[[cell_id_col]]), cell_type_cols, drop = FALSE]
sobj@meta.data <- cbind(sobj@meta.data, celltypes)

log$info(paste0("NOTE: Seurat object identity set to: ", ident))
Idents(sobj) <- ident

log$info(paste0("Saving Seurat object with cell type annotations ..."))
biopipen.utils::save_obj(sobj, outfile)
