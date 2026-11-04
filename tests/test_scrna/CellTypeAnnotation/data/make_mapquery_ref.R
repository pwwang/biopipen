#!/usr/bin/env Rscript
# Build the toy MapQuery reference used by the CellTypeAnnotation tests.
#
# Usage: Rscript make_mapquery_ref.R <outfile>
#
# The reference must provide, for `MapQueryArgs$reference.reduction = "refDR"` and
# `MapQueryArgs$reduction.model = "refUMAP"`:
#   * an SCT assay (default) so the query is normalized with the reference model
#   * `refDR`  : the reference embedding used to find anchors. It MUST carry
#                feature loadings -- ProjectCellEmbeddings takes the anchor
#                features from Loadings(reference[[reference.reduction]]), so a
#                bare CreateDimReducObject() reference fails deep inside
#                FetchResiduals_reference with "incorrect number of dimensions"
#   * `refUMAP`: the model UMAP the query is projected with (needs `return.model`)
#   * a `celltype.l2` label column
# Derived from Seurat's bundled pbmc_small: no download, deterministic.
# The .rds is gitignored (~350 KB binary fixture is not shipped).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: make_mapquery_ref.R <outfile>")
outfile <- args[1]

suppressMessages({
    library(Seurat)
})

set.seed(8525)
ref <- suppressWarnings(UpdateSeuratObject(SeuratObject::pbmc_small))
ref <- NormalizeData(ref, verbose = FALSE)

# two reference labels, mirroring what the tests expect in `celltype.l2`
Idents(ref) <- "letter.idents"
ref$celltype.l2 <- ifelse(Idents(ref) == "A", "B cells", "CD4 T cells")

ref <- SCTransform(ref, verbose = FALSE)

# refDR: PCA kept under its own name so the loadings travel with it
ref <- RunPCA(
    ref,
    npcs = 50,
    reduction.name = "refDR",
    reduction.key = "refdr_",
    verbose = FALSE
)

# refUMAP must be built FROM refDR and carry the learned model, otherwise
# MapQuery cannot project a query onto it
ref <- RunUMAP(
    ref,
    reduction = "refDR",
    dims = 1:30,
    reduction.name = "refUMAP",
    reduction.key = "refUMAP_",
    return.model = TRUE,
    verbose = FALSE
)

# anchors are looked up in VariableFeatures(reference) when reference.reduction
# is used; keep them all explicit
VariableFeatures(ref) <- rownames(ref)

stopifnot(
    !is.null(ref[["refUMAP"]]@misc$model),
    nrow(Loadings(ref[["refDR"]])) > 0,
    length(unique(ref$celltype.l2)) >= 2
)

saveRDS(ref, outfile)
cat(
    "wrote", outfile, ":", paste(dim(ref), collapse = "x"), "cells x features |",
    "reductions:", paste(Reductions(ref), collapse = ","), "|",
    "refDR loadings:", paste(dim(Loadings(ref[["refDR"]])), collapse = "x"), "|",
    "celltype.l2:", paste(unique(ref$celltype.l2), collapse = ","), "\n"
)
