library(loomR)
library(DropletUtils)
library(Matrix)

loomfile <- {{in.loomfile | r}}
outdir <- {{out.outdir | r}}

lfile <- connect(filename = loomfile, mode = "r")

# Extract the expression matrix (genes x cells)
expr_matrix <- t(lfile[["matrix"]][, ])
if (!inherits(expr_matrix, "dgCMatrix")) {
    expr_matrix <- Matrix::Matrix(expr_matrix, sparse = TRUE)
}

# Extract gene names and IDs
gene_names <- lfile[["row_attrs/Gene"]][]

gene_ids <- tryCatch({
   lfile[["row_attrs/GeneID"]][]
}, error = function(e) {
   NULL
})

if (is.null(gene_ids)) {
  gene_ids <- gene_names
}

# Extract cell barcodes
cell_barcodes <- lfile[["col_attrs/CellID"]][]

# Close the LOOM file connection
lfile$close_all()

# Create a data frame for gene information
gene_info <- data.frame(
    gene_id = gene_ids,
    gene_name = gene_names
)

# Write the data to 10X format

write10xCounts(
    path = outdir,
    x = expr_matrix,
    gene.id = gene_info$gene_id,
    gene.symbol = gene_info$gene_name,
    barcodes = cell_barcodes,
    version = "3",
    overwrite = TRUE
)
