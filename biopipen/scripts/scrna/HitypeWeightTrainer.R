# HitypeWeightTrainer.R — train marker weights for hitype from a Seurat
# object via `hitype::train_weights()`. The markers are either given by
# `in.markerfile` (a universal marker table or a native hitype/ScType
# db-format file) or discovered from the data by `hitype::find_markers()`.
# The output is a weighted universal marker table
# (`cell_type`/`gene`/`direction`/`weight`/`level`, one row per marker gene)
# that can be fed to `CellTypeAnnotation` with `tool = "hitype"` /
# `envs.hitype.db` for classification. Requires hitype >= 0.0.6 (universal
# marker format support and the `seed` argument of `train_weights()`).

library(Seurat)
library(hitype)
library(biopipen.utils)

if (packageVersion("hitype") < "0.0.6") {
    stop(paste0(
        "HitypeWeightTrainer requires hitype >= 0.0.6 (the universal marker ",
        "format and the `seed` argument of `train_weights()`). Installed: ",
        as.character(packageVersion("hitype"))
    ))
}

srtobj <- {{in.srtobj | r}}
markers <- {{in.markerfile | r}}
outfile <- {{out.outfile | r}}
ident <- {{envs.ident | r}}
assay <- {{envs.assay | r}}
level <- {{envs.level | r}}
species <- {{envs.species | r}}
cancer <- {{envs.cancer | r}}
tissue <- {{envs.tissue | r}}
find_markers_args <- {{ envs.find_markers | r }}
train_args <- {{ envs.train_weights | r }}

log <- get_logger()

# Source shared marker helpers
biopipen_dir <- {{ biopipen_dir | r }}
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-markers.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-markers.R"))

log$info("Reading the Seurat object ...")
sobj <- read_obj(srtobj)

if (!is.null(ident)) {
    if (!ident %in% colnames(sobj@meta.data)) {
        stop(paste0(
            "The column `", ident, "` (from `envs.ident`) does not exist ",
            "in the meta.data of the Seurat object."
        ))
    }
    log$info("Setting Idents to the `{ident}` column ...")
    Idents(sobj) <- sobj@meta.data[[ident]]
}
if (!is.null(assay)) {
    log$info("Setting the default assay to `{assay}` ...")
    DefaultAssay(sobj) <- assay
}

# The markers, as the `path_to_gs` argument of `train_weights()`: a
# data.frame or a path, in either the universal format (long table with
# `cell_type` and `gene` columns, `weight`/`direction` optional) or the
# native hitype/ScType db format — `gs_prepare()` (hitype >= 0.0.6)
# auto-detects both. The training cell types (marker `cell_type`s, or the
# `cellName`s of a native file) must all exist in the current Idents of the
# Seurat object; `train_weights()` stops otherwise.
if (is.null(markers)) {
    # Discover the markers from the data (per cell type = current Idents)
    if (!is.null(species) || !is.null(cancer) || !is.null(tissue)) {
        stop(paste0(
            "`envs.species`/`envs.cancer`/`envs.tissue` only filter a ",
            "universal marker table given by `in.markerfile`; they are not ",
            "applicable when the markers are discovered from the data by ",
            "`hitype::find_markers()`."
        ))
    }
    unknown_args <- setdiff(names(find_markers_args), formalArgs(find_markers))
    if (length(unknown_args) > 0) {
        stop(paste0(
            "Unknown arguments in `envs.find_markers`: ",
            paste(unknown_args, collapse = ", "),
            " (not arguments of `hitype::find_markers()`)"
        ))
    }
    find_markers_args <- find_markers_args[
        !vapply(find_markers_args, is.null, logical(1))
    ]
    find_args <- c(
        list(exprs = sobj), find_markers_args, list(level = level)
    )
    log$info("Finding markers with `hitype::find_markers()` ...")
    path_to_gs <- do_call(find_markers, find_args)
    counts <- table(path_to_gs$cell_type)
    log$info(
        "Found markers for {length(counts)} cell types: ",
        "{paste(names(counts), as.integer(counts), sep = ' = ', collapse = '; ')}"
    )
} else {
    if (!file.exists(markers)) {
        stop(paste0("The marker file does not exist: ", markers))
    }
    db_markers <- load_marker_table(markers)
    if (is.character(db_markers)) {
        # Native ScType xlsx passthrough (gs_prepare reads it and
        # auto-detects the format)
        stop_on_filtering_native_db(tissue = tissue, cancer = cancer, species = species)
        path_to_gs <- markers
        log$info("Using the native marker file as is: {markers}")
    } else if (is.data.frame(db_markers) && is_marker_canonical(db_markers)) {
        # Universal marker table (a `weight` column is kept; hitype >= 0.0.6
        # reads it)
        path_to_gs <- apply_marker_filters(
            db_markers, tissue = tissue, cancer = cancer, species = species
        )
        log$info("Using the universal marker table: {markers}")
    } else if (is.data.frame(db_markers)) {
        # Native hitype/ScType db-format table (tissueType/cellName/...)
        # without cell_type/gene columns
        stop_on_filtering_native_db(tissue = tissue, cancer = cancer, species = species)
        path_to_gs <- db_markers
        log$info("Using the native db-format marker table: {markers}")
    } else {
        stop(paste0(
            "Cannot recognize the marker format of: ", markers, " ",
            "Use a universal marker table (with `cell_type` and `gene` ",
            "columns) or a native hitype/ScType db-format file."
        ))
    }
}

unknown_args <- setdiff(names(train_args), formalArgs(train_weights))
if (length(unknown_args) > 0) {
    stop(paste0(
        "Unknown arguments in `envs.train_weights`: ",
        paste(unknown_args, collapse = ", "),
        " (not arguments of `hitype::train_weights()`)"
    ))
}
train_args <- train_args[!vapply(train_args, is.null, logical(1))]

log$info("Training marker weights with `hitype::train_weights()` ...")
trained <- do_call(
    train_weights,
    c(list(path_to_gs = path_to_gs, exprs = sobj), train_args, list(level = level))
)

expected_cols <- c("cell_type", "gene", "direction", "weight", "level")
if (!is.data.frame(trained) || !all(expected_cols %in% colnames(trained))) {
    stop(paste0(
        "Unexpected output from `hitype::train_weights()`; expected columns: ",
        paste(expected_cols, collapse = ", ")
    ))
}
log$info(
    "Trained weights for {length(unique(trained$cell_type))} cell types ",
    "({nrow(trained)} markers in total)"
)
if (nrow(trained) > 0) {
    log$info(
        "Weight range: {paste(format(range(trained$weight), digits = 3), collapse = ' - ')}"
    )
}

# Plain write.table, not biopipen's write_table: no `#`-comment header line,
# so the file is re-read with its column names intact by gs_prepare()
# (via read.table) when used as `envs.hitype.db` of `CellTypeAnnotation`.
write.table(trained, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
log$info("Marker weights saved to: {outfile}")
log$info(
    "Use the file with `CellTypeAnnotation` as ",
    "`tool = 'hitype'` / `envs.hitype.db` for classification."
)
