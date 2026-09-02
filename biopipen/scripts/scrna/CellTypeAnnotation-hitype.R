# CellTypeAnnotation-hitype.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_hitype <- function(sobj, ident, tissue, cancer, species, db) {
    library(hitype)

    log <- get_logger()

    if (is.null(db)) { stop("`envs.hitype.db` is not set") }

    # prepare gene sets
    log$info("Preparing gene sets...")
    if (startsWith(db, "hitypedb_") && !grepl(".", db, fixed = TRUE)) {
        # Built-in databases have no cancer/species columns
        stop_on_filtering_native_db(NULL, cancer, species)
        gs_list <- gs_prepare(eval(as.symbol(db)), tissue)
    } else {
        db_markers <- load_marker_table(db)
        if (!is.data.frame(db_markers) || !is_marker_canonical(db_markers)) {
            # Native ScType xlsx/TSV/RDS formats have no cancer/species columns
            # (tissue is still handled via gs_prepare below)
            stop_on_filtering_native_db(NULL, cancer, species)
        }
        if (is.character(db_markers)) {
            # native ScType xlsx passthrough
            gs_list <- gs_prepare(db_markers, tissue)
        } else {
            if (!is.data.frame(db_markers)) {
                stop("Cannot recognize the hitype database format. ",
                     "Use a ScType xlsx/TSV, RDS data.frame, or a universal marker table.")
            }
            if (is_marker_canonical(db_markers)) {
                if (!is.null(tissue) && !"tissue" %in% colnames(db_markers)) {
                    stop(paste0(
                        "`envs.hitype.tissue` is set to `", tissue,
                        "` but the marker table has no `tissue` column."
                    ))
                }
                db_markers <- markers_to_sctype_df(db_markers, tissue, cancer, species)
            }
            # gs_prepare accepts a data.frame directly
            gs_list <- gs_prepare(db_markers, tissue)
        }
    }

    # run RunHitype
    log$info("Running RunHitype...")
    sobj <- RunHitype(sobj, gs_list, threshold = 0.0, make_unique = TRUE)

    log$info("Extracting cell type labels...")
    hitype_labels <- sobj@meta.data %>%
        distinct(!!sym(ident), hitype)
    hitype_labels <- stats::setNames(
        as.list(hitype_labels$hitype),
        hitype_labels[[ident]]
    )

    list(mapping = hitype_labels)
}
