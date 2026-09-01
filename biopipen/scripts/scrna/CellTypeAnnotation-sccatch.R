# CellTypeAnnotation-sccatch.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_sccatch <- function(sobj, ident, sccatch_args) {
    library(scCATCH)

    log <- get_logger()

    if (!is.null(sccatch_args$marker)) {
        cellmatch <- load_marker_table(sccatch_args$marker)
        if (is_marker_canonical(cellmatch)) {
            cellmatch <- markers_to_sccatch_df(cellmatch)
        } else if (!is.data.frame(cellmatch)) {
            stop("The custom marker file for scCATCH must be a table or data.frame.")
        }
        sccatch_args$if_use_custom_marker <- TRUE

        # when if_use_custom_marker is TRUE, sccatch won't filter the markers
        if (!is.null(sccatch_args$species) && "species" %in% colnames(cellmatch)) {
            if (!sccatch_args$species %in% unique(cellmatch$species)) {
                stop(paste(
                    "The species specified in the config does not match the species in the marker file.",
                    "Please check the config file and marker file."
                ))
            }
            cellmatch <- cellmatch[cellmatch$species == sccatch_args$species, , drop = FALSE]
        }

        if (!is.null(sccatch_args$cancer) && "cancer" %in% colnames(cellmatch)) {
            if (!sccatch_args$cancer %in% unique(cellmatch$cancer)) {
                stop(paste(
                    "The cancer type specified in the config does not match the cancer type in the marker file.",
                    "Please check the config file and marker file."
                ))
            }
            cellmatch <- cellmatch[cellmatch$cancer == sccatch_args$cancer, , drop = FALSE]
        }

        if (!is.null(sccatch_args$tissue) && "tissue" %in% colnames(cellmatch)) {
            if (!sccatch_args$tissue %in% unique(cellmatch$tissue)) {
                stop(paste(
                    "The tissue type specified in the config does not match the tissue type in the marker file.",
                    "Please check the config file and marker file."
                ))
            }
            cellmatch <- cellmatch[cellmatch$tissue == sccatch_args$tissue, , drop = FALSE]
        }
    } else {
        sccatch_args$cancer <- sccatch_args$cancer %||% "Normal"
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
