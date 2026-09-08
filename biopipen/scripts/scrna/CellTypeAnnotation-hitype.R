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
                # A universal marker table — consumed natively by hitype
                # (>= 0.0.6). Filter the rows by
                # `envs.tissue`/`envs.cancer`/`envs.species` here and keep
                # the table as is: notably a numeric `weight` column must
                # survive (markers_to_sctype_df(), the sctype route, would
                # drop it).
                if (packageVersion("hitype") < "0.0.6") {
                    stop(paste0(
                        "Universal marker tables with `tool = 'hitype'` ",
                        "require hitype >= 0.0.6 (installed: ",
                        as.character(packageVersion("hitype")), "). ",
                        "Install the latest hitype or use a native ",
                        "db-format file."
                    ))
                }
                db_markers <- apply_marker_filters(
                    db_markers, tissue = tissue, cancer = cancer, species = species
                )
                # Tissues already filtered above; gs_prepare accepts a
                # data.frame directly
                gs_list <- gs_prepare(db_markers, NULL)
            } else {
                # Native hitype/ScType db-format data.frame
                gs_list <- gs_prepare(db_markers, tissue)
            }
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
