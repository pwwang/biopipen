# CellTypeAnnotation-sccatch.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_sccatch <- function(sobj, ident, sccatch_args) {
    library(scCATCH)

    log <- get_logger()

    if (!is.null(sccatch_args$marker)) {
        cellmatch <- read_obj(sccatch_args$marker)
        sccatch_args$if_use_custom_marker <- TRUE
    }
    sccatch_args$marker <- cellmatch

    if (is.integer(sccatch_args$use_method)) {
        sccatch_args$use_method <- as.character(sccatch_args$use_method)
    }

    # Check if there is less than 2 clusters
    num_clusters <- length(unique(sobj@meta.data[[ident]]))
    if (num_clusters < 2) {
        stop(paste(
            "The number of clusters is less than 2.",
            "Sccatch requires at least 2 clusters to perform cell type annotation."
        ))
    }

    log$info("Running createscCATCH ...")
    obj <- createscCATCH(
        data = GetAssayData(sobj, assay = sccatch_args$assay),
        cluster = as.character(sobj@meta.data[[ident]])
    )
    sccatch_args$object <- obj

    log$info("Running findmarkergene ...")
    obj <- do_call(findmarkergene, sccatch_args)

    log$info("Running findcelltype ...")
    obj <- findcelltype(object = obj)

    celltypes <- as.list(obj@celltype$cell_type)
    names(celltypes) <- obj@celltype$cluster

    if (length(celltypes) == 0) {
        log$warn("- No cell types annotated from the database!")
    }

    list(mapping = celltypes)
}
