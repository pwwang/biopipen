# GarnettClassifierTrainer.R — train a `garnett_classifier` from a Seurat
# object and marker genes via `garnett::train_cell_classifier()`.
# The classifier file can be fed to `CellTypeAnnotation` with
# `tool = "garnett"` / `envs.garnett.classifier` for classification.

library(Seurat)
library(monocle3)
library(garnett)
library(SeuratWrappers)
library(biopipen.utils)

srtobj <- {{in.srtobj | r}}
markerfile <- {{in.markerfile | r}}
outfile <- {{out.outfile | r}}
species <- {{envs.species | r}}
cancer <- {{envs.cancer | r}}
tissue <- {{envs.tissue | r}}
db <- {{envs.db | r}}
cds_gene_id_type <- {{envs.cds_gene_id_type | r}}
marker_file_gene_id_type <- {{envs.marker_file_gene_id_type | r}}
classifier_gene_id_type <- {{envs.classifier_gene_id_type | r}}
assay <- {{envs.assay | r}}
min_observations <- {{envs.min_observations | r}}
max_training_samples <- {{envs.max_training_samples | r}}
num_unknown <- {{envs.num_unknown | r}}
propogate_markers <- {{envs.propogate_markers | r}}
lambdas <- {{envs.lambdas | r}}
cores <- {{envs.cores | r}}
seed <- {{envs.seed | r}}
fail_on_degenerate <- {{envs.fail_on_degenerate | r}}

log <- get_logger()
if (!is.null(seed)) { set.seed(seed) } else {
    log$warn(paste(
        "`envs.seed` is not set: garnett drops rare classes and draws the",
        "`Unknown`/training cells at random, so the trained classifier is not",
        "reproducible. Set `envs.seed` for a reproducible run."
    ))
}

# Source shared helper functions
biopipen_dir <- {{ biopipen_dir | r }}
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-markers.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-markers.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-garnett.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-garnett.R"))

# Same glmnet >= 4.0 workaround as the classification side: subtype models
# are assigned via make_predictions() during training too
patch_garnett_make_predictions(log)

# Let garnett's marker-file parser accept Unicode cell type names (e.g.
# "γδ-T cells" from ScTypeDB) instead of failing with "Syntax error 'γ'"
patch_garnett_marker_lexer(log)

# Diagnostics of a trained classifier: one row per fitted node model, with
# `degenerate` TRUE when the selected lambda yields no non-zero coefficient
# (a null model predicting the class priors only). Tolerant of the layout: a
# `garnett_classifier` (igraph tree, fits in the `model` vertex attribute) or
# a plain `list(models = <list of fits>)`. An unexpected layout yields NA
# rows plus a warning instead of an error.
classifier_diagnostics <- function(classifier) {
    na_row <- function(node) data.frame(
        node = node,
        n_classes = NA_integer_,
        classes = NA_character_,
        lambda_min = NA_real_,
        max_abs_coef = NA_real_,
        n_nonzero = NA_integer_,
        degenerate = NA,
        stringsAsFactors = FALSE
    )
    rows <- tryCatch({
        if (isS4(classifier) && is(classifier, "garnett_classifier")) {
            tree <- classifier@classification_tree
            nodes <- igraph::V(tree)$name
            models <- igraph::vertex_attr(tree, "model")
        } else if (is.list(classifier) && !is.null(classifier$models)) {
            models <- classifier$models
            nodes <- names(models) %||% paste0("node", seq_along(models))
        } else {
            stop(
                "unknown classifier layout: ",
                paste(class(classifier), collapse = "/")
            )
        }
        # only the nodes that actually carry a fitted model (in garnett's flat
        # trees that is the root; the other nodes only dispatch on it)
        fitted <- !vapply(models, is.null, logical(1))
        models <- models[fitted]
        nodes <- nodes[fitted]
        lapply(seq_along(models), function(i) {
            model <- models[[i]]
            # cv.glmnet wraps the fitted glmnet model; a plain glmnet fit does
            # not, and then only `lambda.min` (if any) is missing
            fit <- model$glmnet.fit %||% model
            if (is.null(fit$beta)) {
                stop("no `beta` in the model of node ", nodes[i])
            }
            # multinomial fits store one matrix per class, binomial a matrix
            betas <- if (is.list(fit$beta)) fit$beta else list(fit$beta)
            classes <- fit$classnames %||% names(fit$beta)
            lambda_min <- model$lambda.min %||% NA_real_
            if (!is.na(lambda_min) && !is.null(fit$lambda)) {
                j <- which.min(abs(fit$lambda - lambda_min))
                betas <- lapply(betas, function(b) b[, j, drop = FALSE])
            }
            coefs <- unlist(lapply(betas, function(b) as.numeric(as.matrix(b))))
            data.frame(
                node = nodes[i],
                n_classes = length(classes),
                classes = paste(classes, collapse = ","),
                lambda_min = lambda_min,
                max_abs_coef = if (length(coefs)) max(abs(coefs)) else NA_real_,
                n_nonzero = if (length(coefs)) sum(coefs != 0) else NA_integer_,
                degenerate = length(coefs) > 0 && all(coefs == 0),
                stringsAsFactors = FALSE
            )
        })
    }, error = function(e) {
        log$warn(paste0(
            "Could not diagnose the classifier (unexpected layout?): ",
            conditionMessage(e), ". Reporting NA diagnostics."
        ))
        list(na_row("<unknown>"))
    })
    do.call(rbind, rows)
}

# The marker file is either a garnett-native file (passed through unchanged,
# so `subtype of:` hierarchies and `expressed above/below:` rules work) or a
# universal marker table (auto-converted; negative markers become
# `not expressed:` rules). Decided by content (`>` block headers).
if (is_garnett_native_marker(markerfile)) {
    # A native file has no species/cancer/tissue columns to filter by
    stop_on_filtering_native_db(tissue = tissue, cancer = cancer, species = species)
    marker_file <- markerfile
    log$info("Using the garnett-native marker file as is: {basename(markerfile)}")
} else {
    marker_df <- load_marker_table(markerfile)
    if (!is_marker_canonical(marker_df)) {
        stop(paste0(
            "The marker file is neither a garnett-native marker file ",
            "(`> <cell type>` blocks with `expressed:`/`not expressed:` rules, ",
            "see https://cole-trapnell-lab.github.io/garnett/docs/) nor a ",
            "universal marker table (a table with `cell_type` and `gene` ",
            "columns, see the docs of `CellTypeAnnotation`). Got: ", markerfile
        ))
    }
    marker_file <- file.path(
        dirname(outfile),
        paste0(tools::file_path_sans_ext(basename(outfile)), ".markers.txt")
    )
    markers_to_garnett_file(
        marker_df, marker_file,
        tissue = tissue, cancer = cancer, species = species
    )
    log$info("Converted the universal marker table to: {marker_file}")
}

log$info("Reading the Seurat object ...")
sobj <- read_obj(srtobj)

# Gene ID conversion database: "none" or an AnnotationDb/OrgDb package name
db <- db %||% "none"
if (!identical(db, "none")) {
    if (!requireNamespace(db, quietly = TRUE)) {
        stop(paste0(
            "The gene ID database package `", db,
            "` from `envs.db` is not installed."
        ))
    }
    library(db, character.only = TRUE)
    db <- getExportedValue(db, db)
}

log$info("Converting the Seurat object to a monocle3 cell_data_set ...")
# as.cell_data_set is the monocle3 generic; SeuratWrappers only registers
# the method for Seurat objects (loaded above), so call it unqualified
cds <- if (is.null(assay)) {
    as.cell_data_set(sobj)
} else {
    as.cell_data_set(sobj, assay = assay)
}

# garnett's two stages must agree on scale. train_cell_classifier() divides
# counts(cds) by colData(cds)$Size_Factor as supplied, while classify_cells()
# ignores that column and recomputes
#   sfs <- colSums(counts) / (classifier@cell_totals * median(num_genes_expressed))
# i.e. the monocle3 convention (colSums / G, where G = exp(mean(log(colSums)))).
# SeuratWrappers::as.cell_data_set() sets Size_Factor = colSums(counts), which
# satisfies the NULL/NA assert but leaves every linear predictor exactly G-fold
# too large at predict time: exp() overflows, probabilities become NaN, the
# rank_prob_ratio gate then rejects every cell and every cell is silently
# "Unknown" -- with no error, in every arm, however healthy the model is.
# So install monocle3 size factors unconditionally, and verify it took effect.
cds <- monocle3::estimate_size_factors(cds)
.sf <- as.numeric(colData(cds)$Size_Factor)
.ct <- as.numeric(Matrix::colSums(counts(cds)))
.g <- exp(mean(log(.ct[.ct > 0])))
if (anyNA(.sf) || !all(is.finite(.sf))) {
    stop("Could not install finite size factors on the cell_data_set.")
}
if (!isTRUE(all.equal(.g, 1)) &&
    isTRUE(all.equal(.sf, .ct, check.attributes = FALSE))) {
    stop(paste(
        "Garnett requires monocle3 (median-ratio) size factors, but",
        "colData(cds)$Size_Factor still equals colSums(counts) (the",
        "SeuratWrappers::as.cell_data_set() convention). classify_cells()",
        "scores with colSums / (cell_totals * median(num_genes_expressed)),",
        "so training on colSums-scaled values makes the prediction matrix a",
        "factor of G = exp(mean(log(colSums))) too large: probabilities",
        "overflow to NaN and every cell silently becomes 'Unknown'.",
        "Run monocle3::estimate_size_factors(cds) on the converted",
        "cell_data_set before training."
    ))
}

# train_cell_classifier() internally normalizes counts(cds); an assay with
# only a data layer (converted counts are log values) would silently give
# garbage
if (is(counts(cds), "dgCMatrix") && any(counts(cds)@x %% 1 != 0)) {
    log$warn(paste(
        "The counts of the converted cell_data_set are not integers;",
        "garnett::train_cell_classifier() expects raw counts. Results may be",
        "unreliable. Set `envs.assay` to an assay with a counts layer."
    ))
}

# Arguments of train_cell_classifier() (return_initial_assign is deliberately
# not exposed: a classifier must be returned). NULL values use the defaults.
train_args <- list(
    cds_gene_id_type = cds_gene_id_type,
    marker_file_gene_id_type = marker_file_gene_id_type,
    classifier_gene_id_type = classifier_gene_id_type,
    min_observations = min_observations,
    max_training_samples = max_training_samples,
    num_unknown = num_unknown,
    propogate_markers = propogate_markers,
    cores = cores,
    lambdas = lambdas
)
train_args <- train_args[!vapply(train_args, is.null, logical(1))]

log$info("Training the Garnett classifier ...")
classifier <- do_call(
    train_cell_classifier,
    c(list(cds = cds, marker_file = marker_file, db = db), train_args)
)

log$info(
    "Classifier cell types: ",
    "{paste(igraph::V(classifier@classification_tree)$name, collapse = ', ')}"
)

diag <- classifier_diagnostics(classifier)
log$info(
    "Classifier diagnostics (per node model):\n",
    "{paste(utils::capture.output(print(diag, row.names = FALSE)), collapse = '\n')}"
)
diag_file <- paste0(outfile, ".diagnostics.tsv")
utils::write.table(diag, diag_file, sep = "\t", quote = FALSE, row.names = FALSE)
log$info("Classifier diagnostics written to: {diag_file}")

if (nrow(diag) > 0 && isTRUE(all(diag$degenerate))) {
    degenerate_msg <- paste(
        "The fitted models of the classifier have all-zero coefficients: a",
        "null model predicting the class priors only.",
        "Cause: garnett fits every node with glmnet::cv.glmnet() over its",
        "hard-coded lambda grid",
        "`unique(c(1e+05, 50000, seq(10000, 100, by = -200), ...))` with",
        "`standardize = FALSE`, and the cross-validation can select the first",
        "grid point (lambda = 1e+05), where every coefficient is 0.",
        "Every cell will then be classified as `Unknown` by the",
        "`rank_prob_ratio` gate.",
        "Remedy: pass `envs.lambdas`, e.g. `10^seq(1, -3, length.out = 30)`,",
        "so the models are fitted on a sensible lambda grid."
    )
    if (fail_on_degenerate) {
        stop(degenerate_msg)
    } else {
        log$warn(degenerate_msg)
    }
}

save_obj(classifier, outfile)
log$info("Classifier saved to: {outfile}")
