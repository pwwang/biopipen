# CellTypeAnnotation-garnett.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R and GarnettClassifierTrainer.R

# Workaround for garnett 0.2.22 (monocle3 branch) against glmnet >= 4.0:
# `predict()` on a multinomial `cv.glmnet` returns a 3D array, and the
# `as.data.frame()` flattening inside `garnett::make_predictions()` mangles
# the per-class names with a ".lambda.min" suffix (e.g. "T cells.lambda.min").
# The per-type assignment masks then never match the plain cell-type names
# of the classification tree, and every cell comes back "Unknown". Strip
# the suffix from the returned mask names (a no-op once garnett is fixed).
# Idempotent — safe to call before both classify_cells() and
# train_cell_classifier() (subtype models use make_predictions too).
patch_garnett_make_predictions <- function(log) {
    tryCatch({
        ns <- asNamespace("garnett")
        orig_make_predictions <- get("make_predictions", envir = ns)
        monkey_patch("garnett", "make_predictions", function(...) {
            res <- orig_make_predictions(...)
            names(res) <- sub("\\.lambda\\..*$", "", names(res))
            res
        })
    }, error = function(e) {
        log$warn(paste(
            "Failed to patch garnett::make_predictions for glmnet >= 4.0:",
            conditionMessage(e)
        ))
    })
}

annotate_garnett <- function(sobj, ident, garnett_args) {
    library(monocle3)
    library(garnett)
    library(SeuratWrappers)

    log <- get_logger()
    patch_garnett_make_predictions(log)

    classifier_path <- garnett_args$classifier
    if (is.null(classifier_path)) { stop("`envs.garnett.classifier` is not set") }
    if (!file.exists(classifier_path)) {
        stop(paste0("Garnett classifier file does not exist: ", classifier_path))
    }

    log$info("Loading Garnett classifier ...")
    classifier <- read_obj(classifier_path)
    if (!inherits(classifier, "garnett_classifier")) {
        stop(paste0(
            "The file in `envs.garnett.classifier` does not contain a ",
            "Garnett classifier (garnett_classifier) object: ", classifier_path
        ))
    }

    # Gene ID conversion database: "none" or an AnnotationDb/OrgDb package name
    db <- garnett_args$db %||% "none"
    if (!identical(db, "none")) {
        if (!requireNamespace(db, quietly = TRUE)) {
            stop(paste0(
                "The gene ID database package `", db,
                "` from `envs.garnett.db` is not installed."
            ))
        }
        library(db, character.only = TRUE)
        db <- getExportedValue(db, db)
    }

    assay <- garnett_args$assay
    # Keys consumed here must not be forwarded to classify_cells()
    garnett_args$classifier <- NULL
    garnett_args$db <- NULL
    garnett_args$assay <- NULL

    unknown_args <- setdiff(names(garnett_args), formalArgs(classify_cells))
    if (length(unknown_args) > 0) {
        stop(paste0(
            "Unknown arguments in `envs.garnett`: ",
            paste(unknown_args, collapse = ", "),
            " (not arguments of `garnett::classify_cells()`)"
        ))
    }

    log$info("Converting Seurat object to a monocle3 cell_data_set ...")
    # as.cell_data_set is the monocle3 generic; SeuratWrappers only registers
    # the method for Seurat objects (loaded above), so call it unqualified
    cds <- if (is.null(assay)) {
        as.cell_data_set(sobj)
    } else {
        as.cell_data_set(sobj, assay = assay)
    }

    # classify_cells() requires a Size_Factor column (see its source); the
    # value is only a gate, but estimate_size_factors is the canonical way in
    if (is.null(colData(cds)$Size_Factor) || anyNA(colData(cds)$Size_Factor)) {
        log$info("Estimating size factors ...")
        cds <- monocle3::estimate_size_factors(cds)
    }

    # classify_cells() internally normalizes counts(cds); an assay with only a
    # data layer (converted counts are log values) would silently give garbage
    if (is(counts(cds), "dgCMatrix") && any(counts(cds)@x %% 1 != 0)) {
        log$warn(paste(
            "The counts of the converted cell_data_set are not integers;",
            "garnett::classify_cells() expects raw counts. Results may be",
            "unreliable. Set `envs.assay` to an assay with a counts layer."
        ))
    }

    log$info("Classifying cells with Garnett ...")
    classify_args <- garnett_args
    classify_args$cds <- cds
    classify_args$classifier <- classifier
    classify_args$db <- db
    result_cds <- do_call(classify_cells, classify_args)

    labels <- colData(result_cds)$cell_type
    if (is.null(labels)) {
        stop("classify_cells() did not return a `cell_type` column")
    }
    # Cells with zero counts are excluded by classify_cells and get NA
    labels[is.na(labels)] <- "Unknown"

    n_unknown <- sum(labels == "Unknown")
    if (n_unknown == length(labels)) {
        log$warn(paste(
            "All cells are classified as 'Unknown' by Garnett; check",
            "`envs.garnett.db`/`envs.garnett.cds_gene_id_type` (classifier is",
            "trained on", classifier@gene_id_type, "genes)"
        ))
    } else {
        log$info(
            "Garnett classified {length(labels) - n_unknown}/{length(labels)}",
            " cells (Unknown: {n_unknown})"
        )
    }

    result <- data.frame(
        garnett_celltype = unname(labels),
        row.names = colnames(result_cds)
    )

    if (is.null(ident)) {
        list(cell_annotations = result, annotation_col = "garnett_celltype")
    } else {
        log$info("Aggregating Garnett results by cluster ...")
        mapping <- majority_vote(
            labels, as.character(sobj@meta.data[[ident]]),
            unknown = "Unknown"  # garnett's sentinel is capital-U
        )
        list(
            mapping = mapping,
            cell_annotations = result,
            annotation_col = "garnett_celltype"
        )
    }
}
