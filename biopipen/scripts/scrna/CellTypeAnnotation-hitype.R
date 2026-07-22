# CellTypeAnnotation-hitype.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_hitype <- function(sobj, ident, tissue, db) {
    library(hitype)

    log <- get_logger()

    if (is.null(db)) { stop("`hitype_db` is not set") }

    # prepare gene sets
    log$info("Preparing gene sets...")
    if (startsWith(db, "hitypedb_") && !grepl(".", db, fixed = TRUE)) {
        gs_list <- gs_prepare(eval(as.symbol(db)), tissue)
    } else {
        gs_list <- gs_prepare(db, tissue)
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
