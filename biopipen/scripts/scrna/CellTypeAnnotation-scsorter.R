# CellTypeAnnotation-scsorter.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_scsorter <- function(sobj, ident, scsorter_db, scsorter_args) {
    # Check if RunScSorter from the scSorter package is available
    # If not, this is the cran/scSorter package, which doesn't support RunScSorter
    # on Seurat objects
    # We need pwwang/scSorter, which is a fork of scSorter that supports Seurat objects
    if (!requireNamespace("scSorter", quietly = TRUE) ||
        !"RunScSorter" %in% getNamespaceExports("scSorter")) {
        stop(paste(
            "The scSorter package is not installed or does not support RunScSorter.",
            "Please install the pwwang/scSorter package from GitHub."
        ))
    }
    library(scSorter)

    log <- get_logger()

    if (is.null(scsorter_db)) { stop("`envs.scsorter.db` is not set") }

    log$info("Loading scSorter database ...")
    anno <- load_marker_table(scsorter_db)
    if (!is.data.frame(anno)) {
        stop(paste0(
            "scSorter database must be a data.frame, got: ",
            paste(class(anno), collapse = ", ")
        ))
    }
    if (is_marker_canonical(anno)) {
        anno <- markers_to_scsorter_df(
            anno,
            tissue = scsorter_args$tissue,
            cancer = scsorter_args$cancer,
            species = scsorter_args$species
        )
    } else {
        # Native positional tables have no tissue/cancer/species columns
        stop_on_filtering_native_db(
            scsorter_args$tissue, scsorter_args$cancer, scsorter_args$species
        )
        if (ncol(anno) < 2) {
            stop(paste0(
                "scSorter database file must have at least 2 columns: ",
                scsorter_db
            ))
        }
        if (ncol(anno) == 2) {
            colnames(anno) <- c("Type", "Marker")
        } else if (ncol(anno) >= 3) {
            colnames(anno)[1:3] <- c("Type", "Marker", "Weight")
        }
    }
    # The filter envs are consumed by the marker conversion above; never
    # forward them to RunScSorter()
    scsorter_args$tissue <- NULL
    scsorter_args$cancer <- NULL
    scsorter_args$species <- NULL

    log$info("Running RunScSorter...")
    # Set the active identity to the ident column
    Idents(sobj) <- ident
    scsorter_args$db <- NULL  # db is passed separately as scsorter_db
    scsorter_args$object <- sobj
    scsorter_args$anno <- anno
    scsorter_args$mc.cores <- scsorter_args$mc.cores %||% 1L
    sobj <- do_call(RunScSorter, scsorter_args)

    # RunScSorter stores per-cell predictions in the scSorter_celltype column;
    # aggregate to one type per cluster by majority vote
    log$info("Aggregating scSorter results by cluster...")
    mapping <- majority_vote(
        sobj@meta.data$scSorter_celltype,
        as.character(sobj@meta.data[[ident]])
    )

    list(mapping = mapping)
}
