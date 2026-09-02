library(rlang)
library(dplyr)
library(Seurat)
library(tidyseurat)
library(biopipen.utils)
set.seed(8525)

sobjfile <- {{in.sobjfile | r}}
outfile <- {{out.outfile | r}}
tool <- {{envs.tool | r}}
ident <- {{envs.ident | r }}
merge <- {{envs.merge | r}}
add_prefix <- {{envs.add_prefix | r}}
outtype <- {{envs.outtype | r}}
cases <- {{envs.cases | r}}
ncores <- {{envs.ncores | r}}
anno_col <- {{envs.anno_col | r}}
set_ident <- {{envs.set_ident | r}}

assay <- {{envs.assay | r}}
# Tool-specific envs (new style)
sctype <- {{envs.sctype | r}}
hitype <- {{envs.hitype | r}}
scsorter <- {{envs.scsorter | r}}
scina <- {{envs.scina | r}}
singler <- {{envs.singler | r}}
schdeepinsight <- {{envs.schdeepinsight | r}}
gptcelltype <- {{envs.gptcelltype | r}}
cellassign <- {{envs.cellassign | r}}
scbert <- {{envs.scbert | r}}
cellid <- {{envs.cellid | r}}
sccatch <- {{envs.sccatch | r}}
celltypist <- {{envs.celltypist | r}}

# Deprecated envs (still supported with a warning)
sctype_tissue <- {{ envs["sctype_tissue"] | default: None | r }}
sctype_db <- {{ envs["sctype_db"] | default: None | r }}
hitype_tissue <- {{ envs["hitype_tissue"] | default: None | r }}
hitype_db <- {{ envs["hitype_db"] | default: None | r }}
scsorter_db <- {{ envs["scsorter_db"] | default: None | r }}
scsorter_args <- {{ envs["scsorter_args"] | default: None | r }}
scina_db <- {{ envs["scina_db"] | default: None | r }}
scina_args <- {{ envs["scina_args"] | default: None | r }}
singler_db <- {{ envs["singler_db"] | default: None | r }}
singler_args <- {{ envs["singler_args"] | default: None | r }}
schdeepinsight_ref <- {{ envs["schdeepinsight_ref"] | default: None | r }}
schdeepinsight_args <- {{ envs["schdeepinsight_args"] | default: None | r }}
gptcelltype_args <- {{ envs["gptcelltype_args"] | default: None | r }}
cellassign_db <- {{ envs["cellassign_db"] | default: None | r }}
cellassign_args <- {{ envs["cellassign_args"] | default: None | r }}
scbert_ref <- {{ envs["scbert_ref"] | default: None | r }}
scbert_model <- {{ envs["scbert_model"] | default: None | r }}
scbert_label_dict <- {{ envs["scbert_label_dict"] | default: None | r }}
scbert_args <- {{ envs["scbert_args"] | default: None | r }}
cellid_db <- {{ envs["cellid_db"] | default: None | r }}
cellid_args <- {{ envs["cellid_args"] | default: None | r }}
sccatch_args <- {{ envs["sccatch_args"] | default: None | r }}
celltypist_args <- {{ envs["celltypist_args"] | default: None | r }}
newcol <- {{ envs["newcol"] | default: None | r }}
backup_col <- {{ envs["backup_col"] | default: None | r }}
cell_types <- {{envs.cell_types | r}}
more_cell_types <- {{envs.more_cell_types | r}}

qs2::qopt("nthreads", value = ncores)

log <- get_logger()

# Source all tool function definitions
biopipen_dir <- {{ biopipen_dir | r }}
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-hitype.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-hitype.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "sctype.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "sctype.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-sctype.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-sctype.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-sccatch.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-sccatch.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-celltypist.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-celltypist.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-scsorter.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-scsorter.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-direct.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-direct.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-cell.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-cell.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-scina.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-scina.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-singler.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-singler.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-schdeepinsight.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-schdeepinsight.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-gptcelltype.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-gptcelltype.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-cellassign.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-cellassign.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-scbert.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-scbert.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-cellid.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-cellid.R"))
# {{ biopipen_dir | joinpaths: "scripts", "scrna", "CellTypeAnnotation-markers.R" | getmtime | int }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-markers.R"))

# Majority-vote helper for cell-level tools with ident
majority_vote <- function(labels, clusters, unknown = "unknown") {
    mapping <- list()
    for (cl in unique(as.character(clusters))) {
        cl_labels <- labels[as.character(clusters) == cl]
        cl_labels <- cl_labels[!is.na(cl_labels)]
        known <- cl_labels[cl_labels != unknown]
        mapping[[cl]] <- if (length(known) > 0) {
            names(which.max(table(known)))
        } else {
            unknown
        }
    }
    mapping
}

# Deprecated envs compatibility
DEPRECATED_ENVS <- list(
    sctype = list(tissue = "sctype_tissue", db = "sctype_db"),
    hitype = list(tissue = "hitype_tissue", db = "hitype_db"),
    scsorter = list(db = "scsorter_db", args = "scsorter_args"),
    scina = list(db = "scina_db", args = "scina_args"),
    singler = list(db = "singler_db", args = "singler_args"),
    schdeepinsight = list(
        ref = "schdeepinsight_ref",
        args = "schdeepinsight_args"
    ),
    gptcelltype = list(args = "gptcelltype_args"),
    cellassign = list(db = "cellassign_db", args = "cellassign_args"),
    scbert = list(
        ref = "scbert_ref",
        model = "scbert_model",
        label_dict = "scbert_label_dict",
        args = "scbert_args"
    ),
    cellid = list(db = "cellid_db", args = "cellid_args"),
    sccatch = list(args = "sccatch_args"),
    celltypist = list(args = "celltypist_args")
)

.deprecation_warned <- character(0)

warn_deprecated_once <- function(env_name, message) {
    if (env_name %in% .deprecation_warned) return(invisible(NULL))
    .deprecation_warned <<- c(.deprecation_warned, env_name)
    log$warn(message)
}

# Old-style envs take precedence over new-style ones
normalize_deprecated <- function(case) {
    if (!is.null(case$newcol)) {
        warn_deprecated_once(
            "newcol",
            "`envs.newcol` is deprecated, use `envs.anno_col` instead"
        )
        case$anno_col <- case$newcol
        case$newcol <- NULL
    }
    if (!is.null(case$backup_col)) {
        warn_deprecated_once(
            "backup_col",
            paste0(
                "`envs.backup_col` is deprecated and ignored; ",
                "the original identity column is never modified"
            )
        )
        case$backup_col <- NULL
    }

    for (tool_name in names(DEPRECATED_ENVS)) {
        for (new_key in names(DEPRECATED_ENVS[[tool_name]])) {
            old_key <- DEPRECATED_ENVS[[tool_name]][[new_key]]
            old_val <- case[[old_key]]
            if (is.null(old_val) || length(old_val) == 0) next

            new_ref <- if (identical(new_key, "args")) {
                paste0("envs.", tool_name)
            } else {
                paste0("envs.", tool_name, ".", new_key)
            }
            warn_deprecated_once(
                old_key,
                paste0(
                    "`envs.", old_key, "` is deprecated, use `", new_ref,
                    "` instead"
                )
            )

            ns <- case[[tool_name]] %||% list()
            if (identical(new_key, "args")) {
                ns <- list_update(ns, old_val)
            } else {
                ns[[new_key]] <- old_val
            }
            case[[tool_name]] <- ns
            case[[old_key]] <- NULL
        }
    }

    case
}

# Build defaults list for expand_cases (exclude ncores — not inherited)
defaults <- list(
    tool = tool,
    ident = ident,
    anno_col = anno_col,
    set_ident = set_ident,
    merge = merge,
    sctype = sctype,
    hitype = hitype,
    scsorter = scsorter,
    scina = scina,
    singler = singler,
    schdeepinsight = schdeepinsight,
    gptcelltype = gptcelltype,
    cellassign = cellassign,
    scbert = scbert,
    cellid = cellid,
    sccatch = sccatch,
    celltypist = celltypist,
    cell_types = cell_types,
    more_cell_types = more_cell_types,
    # Deprecated envs (ignored when NULL)
    sctype_tissue = sctype_tissue,
    sctype_db = sctype_db,
    hitype_tissue = hitype_tissue,
    hitype_db = hitype_db,
    scsorter_db = scsorter_db,
    scsorter_args = scsorter_args,
    scina_db = scina_db,
    scina_args = scina_args,
    singler_db = singler_db,
    singler_args = singler_args,
    schdeepinsight_ref = schdeepinsight_ref,
    schdeepinsight_args = schdeepinsight_args,
    gptcelltype_args = gptcelltype_args,
    cellassign_db = cellassign_db,
    cellassign_args = cellassign_args,
    scbert_ref = scbert_ref,
    scbert_model = scbert_model,
    scbert_label_dict = scbert_label_dict,
    scbert_args = scbert_args,
    cellid_db = cellid_db,
    cellid_args = cellid_args,
    sccatch_args = sccatch_args,
    celltypist_args = celltypist_args,
    newcol = newcol,
    backup_col = backup_col
)

cases <- expand_cases(cases, defaults, default_case = "DEFAULT")
cases <- lapply(cases, normalize_deprecated)

# Cluster-based tools
CLUSTER_LEVEL_TOOLS <- c("hitype", "sctype", "sccatch", "singler", "scsorter", "gptcelltype", "direct")
CELL_LEVEL_TOOLS <- c("scina", "cellassign", "cellid", "scbert", "schdeepinsight", "cell", "celltypist")
PYTHON_TOOLS <- c("celltypist", "schdeepinsight", "scbert")

# Handle the edge case: single DEFAULT case with direct tool and empty cell_types
# Backward compat: create a symlink instead of processing
if (length(cases) == 1 &&
    identical(names(cases)[1], "DEFAULT") &&
    identical(cases[["DEFAULT"]]$tool, "direct") &&
    (is.null(cases[["DEFAULT"]]$cell_types) ||
     length(cases[["DEFAULT"]]$cell_types) == 0)) {
    log$warn("No cell types are given and no cases configured!")
    log$info("Creating symbolic link to input file...")
    file.symlink(sobjfile, outfile)
    quit(save = "no", status = 0)
}

log$info("Reading Seurat object...")
sobj <- read_obj(sobjfile)
# Capture the original identity column; RenameSeuratIdents sets Idents() when
# the ident column is the active identity, so restore it after each case and
# let `set_ident` decide the final identity
original_ident <- GetIdentityColumn(sobj)

outdir <- dirname(outfile)
outprefix <- file.path(outdir, tools::file_path_sans_ext(basename(outfile)))

# Pre-convert to h5ad if any case needs it
needs_h5ad <- any(sapply(cases, function(c)
    c$tool %in% PYTHON_TOOLS
))
h5ad_path <- NULL
if (needs_h5ad) {
    if (endsWith(tolower(sobjfile), ".h5ad")) {
        h5ad_path <- sobjfile
    } else {
        log$info("Pre-converting Seurat object to h5ad for Python-based tools...")
        h5ad_path <- paste0(outprefix, ".shared.h5ad")
        ConvertSeuratToAnnData(
            sobj,
            outfile = h5ad_path,
            assay = assay,
            log = log
        )
    }
}

# Resolve ident for each case and determine tool type
cases <- lapply(cases, function(case) {
    case$tool <- tolower(case$tool)
    if (identical(case$ident, "ident")) {
        # "ident" is an alias for the identity column
        case$ident <- GetIdentityColumn(sobj)
    }
    has_ident <- !is.null(case$ident)
    is_cluster_based <- case$tool %in% CLUSTER_LEVEL_TOOLS ||
        (case$tool %in% CELL_LEVEL_TOOLS && has_ident) ||
        (identical(case$tool, "celltypist") &&
            !is.null(case$celltypist$over_clustering) &&
            !isFALSE(case$celltypist$over_clustering))
    if (is_cluster_based && is.null(case$ident)) {
        # Only resolve the identity column for cluster-based cases
        case$ident <- GetIdentityColumn(sobj)
    }
    case$is_cluster_based <- is_cluster_based
    case
})

# Run each case
run_case <- function(case_name) {
    case <- cases[[case_name]]

    log$info("Running annotation for case: {case_name}")
    tool_name <- case$tool
    tool_cfg <- case[[tool_name]] %||% list()
    if (identical(tool_name, "celltypist") && is.null(tool_cfg$assay)) {
        tool_cfg$assay <- assay
    }
    if (identical(tool_name, "scsorter") && is.null(tool_cfg$assay)) {
        tool_cfg$assay <- assay
    }

    result <- switch(tool_name,
        hitype = annotate_hitype(
            sobj, case$ident, tool_cfg$tissue, tool_cfg$cancer, tool_cfg$species, tool_cfg$db
        ),
        sctype = annotate_sctype(
            sobj, case$ident, tool_cfg$tissue, tool_cfg$cancer, tool_cfg$species, tool_cfg$db
        ),
        sccatch = annotate_sccatch(
            sobj, case$ident, tool_cfg
        ),
        celltypist = annotate_celltypist(
            sobj, case$ident, tool_cfg, outdir,
            h5ad_path = h5ad_path, case_id = case_name
        ),
        scsorter = annotate_scsorter(
            sobj, case$ident, tool_cfg$db, tool_cfg
        ),
        scina = annotate_scina(
            sobj, case$ident, tool_cfg$db, tool_cfg
        ),
        singler = annotate_singler(
            sobj, case$ident, tool_cfg$db, tool_cfg
        ),
        schdeepinsight = annotate_schdeepinsight(
            sobj, case$ident, tool_cfg$ref, tool_cfg, outdir,
            h5ad_path = h5ad_path, case_id = case_name
        ),
        gptcelltype = annotate_gptcelltype(
            sobj, case$ident, tool_cfg
        ),
        cellassign = annotate_cellassign(
            sobj, case$ident, tool_cfg$db, tool_cfg
        ),
        scbert = annotate_scbert(
            sobj, case$ident, tool_cfg$ref, tool_cfg$model,
            tool_cfg$label_dict, tool_cfg, outdir,
            h5ad_path = h5ad_path, case_id = case_name
        ),
        cellid = annotate_cellid(
            sobj, case$ident, tool_cfg$db, tool_cfg
        ),
        direct = annotate_direct(
            sobj, case$ident, case$cell_types, case$more_cell_types
        ),
        cell = annotate_cell(sobj, case$cell_types, case$ident),
        stop(paste0("Unknown tool: ", case$tool))
    )

    # For celltypist, the cluster mapping is keyed on the over_clustering column
    result_ident <- case$ident
    if (identical(tool_name, "celltypist") &&
        !is.null(tool_cfg$over_clustering) &&
        !isFALSE(tool_cfg$over_clustering)) {
        result_ident <- tool_cfg$over_clustering
    }

    list(
        name = case_name,
        mapping = result$mapping,
        cell_annotations = result$cell_annotations,
        annotation_col = result$annotation_col,
        more = result$more,
        ident = result_ident,
        anno_col = case$anno_col,
        set_ident = isTRUE(case$set_ident),
        merge = case$merge,
        is_cluster_based = case$is_cluster_based
    )
}

if (ncores > 1 && length(cases) > 1) {
    library(parallel)
    log$info("Running {length(cases)} case(s) in parallel with {ncores} core(s)...")
    results <- mclapply(
        names(cases),
        run_case,
        mc.cores = min(ncores, length(cases)),
        mc.preschedule = FALSE
    )
} else {
    log$info("Running {length(cases)} case(s) sequentially...")
    results <- lapply(names(cases), run_case)
}

# Apply results to Seurat object in order
case_anno_cols <- list()
set_ident_cases <- list()
for (res in results) {
    if (is.null(res)) next

    # Determine prefix for non-DEFAULT cases
    prefix <- ""
    if (!identical(res$name, "DEFAULT") &&
        (is.null(add_prefix) || isTRUE(add_prefix))) {
        prefix <- paste0(res$name, "_")
    }

    # Add per-cell annotations to metadata (for cell-level tools)
    if (!is.null(res$cell_annotations)) {
        cell_annotations <- res$cell_annotations
        # Prefix column names for non-DEFAULT cases
        if (nzchar(prefix)) {
            colnames(cell_annotations) <- paste0(
                prefix, colnames(cell_annotations)
            )
        }

        # Remove existing columns that will be overwritten
        existing <- intersect(
            colnames(cell_annotations),
            colnames(sobj@meta.data)
        )
        if (length(existing) > 0) {
            sobj@meta.data[, existing] <- NULL
        }
        sobj@meta.data <- cbind(sobj@meta.data, cell_annotations)
        log$info("Adding cell-level annotations to metadata")
    }

    if (res$is_cluster_based) {
        case_anno_col <- paste0(prefix, res$anno_col)
        log$info("Adding annotation as new column: {case_anno_col}")
        sobj <- RenameSeuratIdents(
            sobj,
            mapping = res$mapping,
            ident = res$ident,
            save_as = case_anno_col,
            merge = res$merge
        )
        if (!is.null(original_ident)) {
            Idents(sobj) <- original_ident
        }
        case_anno_cols[[res$name]] <- case_anno_col
    } else {
        # Cell-only mode: rename the per-cell annotation column to anno_col
        annotation_col <- paste0(prefix, colnames(res$cell_annotations)[1])
        case_anno_col <- paste0(prefix, res$anno_col)
        if (!identical(annotation_col, case_anno_col)) {
            colnames(sobj@meta.data)[
                colnames(sobj@meta.data) == annotation_col
            ] <- case_anno_col
            log$info(
                "Renamed annotation column '{annotation_col}' to '{case_anno_col}'"
            )
        }
        case_anno_cols[[res$name]] <- case_anno_col
    }

    if (isTRUE(res$set_ident)) {
        set_ident_cases[[length(set_ident_cases) + 1]] <- res$name
    }

    # Handle additional columns from more_cell_types (direct tool)
    if (!is.null(res$more)) {
        for (key in names(res$more)) {
            more_col <- paste0(prefix, key)
            log$info("Adding additional annotation column: {more_col}")
            sobj <- RenameSeuratIdents(
                sobj,
                mapping = res$more[[key]],
                ident = res$ident,
                save_as = more_col,
                merge = res$merge
            )
            if (!is.null(original_ident)) {
                Idents(sobj) <- original_ident
            }
        }
    }
}

# Save cluster-to-cell-type mappings from all cluster-based cases
# Collect mappings per case
case_mappings <- list()
all_clusters <- character(0)
for (res in results) {
    if (is.null(res)) next
    if (res$is_cluster_based && !is.null(res$mapping)) {
        case_mappings[[res$name]] <- res$mapping
        all_clusters <- union(all_clusters, names(res$mapping))
    }
}
if (length(case_mappings) > 0) {
    all_clusters <- sort(all_clusters)
    cluster_df <- data.frame(
        Cluster = all_clusters,
        stringsAsFactors = FALSE
    )
    # Add the size (number of cells) of each cluster
    cluster_sizes <- list()
    for (res in results) {
        if (is.null(res) || !res$is_cluster_based || is.null(res$mapping)) next
        if (!res$ident %in% colnames(sobj@meta.data)) next
        cluster_sizes[[res$name]] <- table(as.character(sobj@meta.data[[res$ident]]))
    }
    if (length(cluster_sizes) > 0) {
        cluster_df$Size <- sapply(all_clusters, function(cl) {
            for (sz in cluster_sizes) {
                if (cl %in% names(sz)) return(as.integer(sz[[cl]]))
            }
            NA_integer_
        })
    }
    for (case_name in names(case_mappings)) {
        cluster_df[[case_name]] <- sapply(all_clusters, function(cl) {
            case_mappings[[case_name]][[cl]] %||% NA_character_
        })
    }
    tsv_file <- paste0(outprefix, ".cluster2celltype.tsv")
    write_table(
        cluster_df, tsv_file,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    log$info("Saved cluster-to-cell-type mappings to: {basename(tsv_file)}")
}

# Save per-cell annotations from all cell-level cases
cell_case_cols <- list()
for (res in results) {
    if (is.null(res) || is.null(res$cell_annotations)) next

    prefix <- ""
    if (!identical(res$name, "DEFAULT") &&
        (is.null(add_prefix) || isTRUE(add_prefix))) {
        prefix <- paste0(res$name, "_")
    }

    if (res$is_cluster_based) {
        # Both-mode: the tool-specific per-cell column is kept
        cell_case_cols[[res$name]] <- paste0(
            prefix, colnames(res$cell_annotations)[1]
        )
    } else {
        # Cell-only mode: the annotation column was renamed to anno_col
        cell_case_cols[[res$name]] <- paste0(prefix, res$anno_col)
    }
}
if (length(cell_case_cols) > 0) {
    cell_df <- data.frame(
        Cell = rownames(sobj@meta.data),
        stringsAsFactors = FALSE
    )
    for (case_name in names(cell_case_cols)) {
        col_name <- cell_case_cols[[case_name]]
        if (col_name %in% colnames(sobj@meta.data)) {
            cell_df[[case_name]] <- sobj@meta.data[[col_name]]
        } else {
            cell_df[[case_name]] <- NA_character_
        }
    }
    tsv_file <- paste0(outprefix, ".cell2celltype.tsv")
    write_table(
        cell_df, tsv_file,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    log$info("Saved per-cell annotations to: {basename(tsv_file)}")
}

# Set identity to the annotation column if set_ident is enabled
if (length(set_ident_cases) > 0) {
    if (length(set_ident_cases) > 1) {
        last_case_name <- set_ident_cases[[length(set_ident_cases)]]
        log$warn(paste0(
            "Multiple cases have set_ident enabled; ",
            "the LAST case (", last_case_name, ") wins"
        ))
    }
    last_case <- set_ident_cases[[length(set_ident_cases)]]
    last_col <- case_anno_cols[[last_case]]
    log$info(paste0(
        "SETTING IDENTS to '", last_col, "' (from case '", last_case, "') ..."
    ))
    Idents(sobj) <- last_col
}

log$info("Saving Seurat object...")
save_obj(sobj, outfile)
