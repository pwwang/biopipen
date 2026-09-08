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

log <- get_logger()
if (!is.null(seed)) { set.seed(seed) }

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

# train_cell_classifier() requires a Size_Factor column (hard assert); the
# value is only a gate, but estimate_size_factors is the canonical way in
if (is.null(colData(cds)$Size_Factor) || anyNA(colData(cds)$Size_Factor)) {
    log$info("Estimating size factors ...")
    cds <- monocle3::estimate_size_factors(cds)
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

save_obj(classifier, outfile)
log$info("Classifier saved to: {outfile}")
