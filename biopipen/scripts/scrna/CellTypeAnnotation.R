library(rlang)
library(dplyr)
library(Seurat)
library(tidyseurat)
library(biopipen.utils)
set.seed(8525)

sobjfile <- {{in.sobjfile | r}}
outfile <- {{out.outfile | r}}
tool <- {{envs.tool | downcase | r}}
ident <- {{envs.ident | r }}
merge <- {{envs.merge | r}}
add_prefix <- {{envs.add_prefix | r}}
outtype <- {{envs.outtype | r}}
cases <- {{envs.cases | r}}
ncores <- {{envs.ncores | r}}
anno_col <- {{envs.anno_col | r}}
set_ident <- {{envs.set_ident | r}}
subset <- {{envs.subset | r}}
post_mutaters <- {{envs.post_mutaters | r}}
report <- {{envs.report | r}}
joboutdir <- {{job.outdir | r}}

assay <- {{envs.assay | r}}
layer <- {{envs.layer | r}}
# Tool-specific envs (new style)
sctype <- {{envs.sctype | r}}
hitype <- {{envs.hitype | r}}
scsorter <- {{envs.scsorter | r}}
scina <- {{envs.scina | r}}
singler <- {{envs.singler | r}}
garnett <- {{envs.garnett | r}}
schdeepinsight <- {{envs.schdeepinsight | r}}
llmcelltype <- {{envs.llmcelltype | r}}
cellassign <- {{envs.cellassign | r}}
scbert <- {{envs.scbert | r}}
cellid <- {{envs.cellid | r}}
sccatch <- {{envs.sccatch | r}}
celltypist <- {{envs.celltypist | r}}
scagenttype <- {{envs.scagenttype | r}}
ucell <- {{envs.ucell | r}}
aucell <- {{envs.aucell | r}}
gsva <- {{envs.gsva | r}}
singscore <- {{envs.singscore | r}}
scmap <- {{envs.scmap | r}}
cheetah <- {{envs.cheetah | r}}
scclassify <- {{envs.scclassify | r}}
scpred <- {{envs.scpred | r}}
azimuth <- {{envs.azimuth | r}}
scsa <- {{envs.scsa | r}}
maca <- {{envs.maca | r}}
scmapnet <- {{envs.scmapnet | r}}
mllmcelltype <- {{envs.mllmcelltype | r}}
lict <- {{envs.lict | r}}
mapquery <- {{envs.mapquery | r}}
scinfer <- {{envs.scinfer | r}}

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
llmcelltype_args <- {{ envs["llmcelltype_args"] | default: None | r }}
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

# Cache directory for the tool conversions (h5ad) and the per-tool scratch files.
# It is per-process on purpose: forked `mclapply` cases share it, so a python
# tool converts the object once per run instead of once per case.
cache_dir <- tempdir()

# The annotation tools live in biopipen.utils::RunCellTypeAnnotation() now; this
# script keeps the process surface: cases, prefixes, applying the mappings to the
# object, the TSV outputs and saving.
biopipen_dir <- {{ biopipen_dir | r }}

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
    llmcelltype = list(args = "llmcelltype_args"),
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
    subset = subset,
    merge = merge,
    sctype = sctype,
    hitype = hitype,
    scsorter = scsorter,
    scina = scina,
    singler = singler,
    garnett = garnett,
    schdeepinsight = schdeepinsight,
    llmcelltype = llmcelltype,
    cellassign = cellassign,
    scbert = scbert,
    cellid = cellid,
    sccatch = sccatch,
    celltypist = celltypist,
    scagenttype = scagenttype,
    ucell = ucell,
    aucell = aucell,
    gsva = gsva,
    singscore = singscore,
    scmap = scmap,
    cheetah = cheetah,
    scclassify = scclassify,
    scpred = scpred,
    azimuth = azimuth,
    scsa = scsa,
    maca = maca,
    scmapnet = scmapnet,
    mllmcelltype = mllmcelltype,
    lict = lict,
    mapquery = mapquery,
    scinfer = scinfer,
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
    llmcelltype_args = llmcelltype_args,
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

# The identity of the input object; the purity of a cell-level annotation is
# computed against it (`seurat_clusters` from SeuratClustering, for example)
default_ident <- GetIdentityColumn(sobj)

outdir <- dirname(outfile)
outprefix <- file.path(outdir, tools::file_path_sans_ext(basename(outfile)))

# Resolve ident for each case and read the tool level from the package's tool table
cases <- lapply(cases, function(case) {
    case$tool <- tolower(case$tool)
    if (identical(case$ident, "ident")) {
        # "ident" is an alias for the identity column
        case$ident <- GetIdentityColumn(sobj)
    }
    desc <- tryCatch(
        celltype_annotation_tools(case$tool)[[1]],
        error = function(e) stop(paste0("Unknown tool: ", case$tool))
    )
    is_cluster_based <- identical(desc$level, "cluster") ||
        !is.null(case$ident) ||
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

# With several mapquery cases, each case's reference projection would be
# written under the same name (ref.umap by default); keep them all by giving
# each case its own reduction name, unless the case set one explicitly
mapquery_cases <- names(cases)[vapply(
    cases, function(case) identical(case$tool, "mapquery"), logical(1)
)]
if (length(mapquery_cases) > 1) {
    for (case_name in mapquery_cases) {
        existing <- cases[[case_name]]$mapquery$map_query$projectumap.args$reduction.name
        cases[[case_name]]$mapquery$map_query$projectumap.args$reduction.name <-
            existing %||% paste0(case_name, ".ref.umap")
    }
    log$info(paste0(
        "Multiple mapquery cases: keeping per-case projections as ",
        paste0("'", mapquery_cases, ".ref.umap'", collapse = ", ")
    ))
}

# Run each case
run_case <- function(case_name) {
    case <- cases[[case_name]]

    log$info("Running annotation for case: {case_name}")
    tool_name <- case$tool
    tool_cfg <- if (tool_name %in% c("direct", "cell")) {
        # `cell_types`/`more_cell_types` are case-level envs, not envs.<tool>
        list(cell_types = case$cell_types, more_cell_types = case$more_cell_types)
    } else {
        case[[tool_name]] %||% list()
    }
    if (tool_name %in% c("celltypist", "garnett", "cellassign")) {
        tool_cfg$assay <- tool_cfg$assay %||% assay
    }
    if (tool_name %in% c(
        "scsorter", "sctype", "hitype", "singler", "scina", "sccatch", "aucell",
        "gsva", "singscore", "scmap", "cheetah", "scclassify"
    )) {
        tool_cfg$assay <- tool_cfg$assay %||% assay
        tool_cfg$layer <- tool_cfg$layer %||% layer
    }

    # Restrict the tool to the case's subset; the cells outside it get NA
    # when the record is embedded
    sobj_case <- sobj
    subset_cells <- NULL
    if (!is.null(case$subset) && length(case$subset) > 0 &&
        any(nzchar(case$subset))) {
        subobj <- filter(sobj, !!parse_expr(case$subset))
        if (ncol(subobj) == 0) {
            stop(paste0(
                "No cells left after applying `subset` for case '", case_name,
                "': ", case$subset
            ))
        }
        log$info(paste0(
            "Subsetting case '", case_name, "' to ", ncol(subobj),
            " of ", ncol(sobj), " cells"
        ))
        subset_cells <- colnames(subobj)
        sobj_case <- subobj
    }

    if (identical(tool_name, "mapquery")) {
        # The mapped object carries the reference projection (ref.umap by
        # default); hand back only its compact payload -- the whole object
        # dies with the fork
        existing <- Reductions(sobj_case)
        res <- RunCellTypeAnnotation(
            sobj_case, tool_name,
            args = tool_cfg, ident = case$ident, cache = cache_dir, log = log,
            return_object = TRUE
        )
        payload <- SeuratProjectionPayload(res$object, exclude = existing)
        record <- res$record
        reductions <- payload
    } else {
        record <- RunCellTypeAnnotation(
            sobj_case, tool_name,
            args = tool_cfg, ident = case$ident, cache = cache_dir, log = log
        )
        reductions <- list()
    }

    # For celltypist, the cluster mapping is keyed on the over_clustering column
    result_ident <- case$ident
    if (identical(tool_name, "celltypist") &&
        !is.null(tool_cfg$over_clustering) &&
        !isFALSE(tool_cfg$over_clustering)) {
        result_ident <- tool_cfg$over_clustering
    }

    # The per-cell scores the tool returned (`record$scores`), when it reports
    # one: a data.frame of numeric columns with cell row names. Anything else
    # is ignored with a warning, so a lagging runner cannot break the run.
    scores <- record$scores
    if (!is.null(scores)) {
        valid_scores <- is.data.frame(scores) &&
            nrow(scores) > 0 &&
            ncol(scores) > 0 &&
            .row_names_info(scores) > 0 &&
            all(vapply(scores, is.numeric, logical(1)))
        if (!valid_scores) {
            log$warn(paste0(
                "Ignoring `scores` of case '", case_name, "': expected a ",
                "data.frame of numeric columns with cell row names"
            ))
            scores <- NULL
        }
    }

    list(
        name = case_name,
        record = record,
        mapping = record$mapping,
        type = record$type,
        cells = record$cells,
        ident = result_ident,
        anno_col = case$anno_col,
        set_ident = isTRUE(case$set_ident),
        merge = case$merge,
        subset_cells = subset_cells,
        is_cluster_based = !is.null(record$mapping),
        reductions = reductions,
        scores = scores
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

# forked cases swallow their error into a "try-error" string; surface it here,
# naming the case, instead of failing later on a cryptic `$` on an atomic vector
failed <- vapply(results, inherits, logical(1), what = "try-error")
if (any(failed)) {
    stop(paste0(
        "[CellTypeAnnotation] case(s) failed: ",
        paste0(
            names(cases)[failed], ": ",
            vapply(
                results[failed],
                function(x) trimws(as.character(x)),
                character(1)
            ),
            collapse = "; "
        )
    ), call. = FALSE)
}

### --- scores: start ---
# Carry a case's per-cell tool scores (`record$scores`) into the metadata,
# mirroring its label columns: the same prefix, and NA for the cells without a
# score or outside the case's subset. `scores` must be a data.frame of numeric
# columns with cell row names; anything else is ignored with a single warning,
# and an existing metadata column is never overwritten.
add_case_scores <- function(sobj, scores, prefix, subset_cells, log) {
    if (is.null(scores) || NCOL(scores) == 0 || NROW(scores) == 0) {
        return(sobj)
    }
    if (!is.data.frame(scores) ||
            .row_names_info(scores) <= 0 ||
            !all(vapply(scores, is.numeric, logical(1)))) {
        log$warn(paste0(
            "Ignoring the scores: expected a data.frame of numeric columns ",
            "with cell row names"
        ))
        return(sobj)
    }
    cells <- colnames(sobj)
    for (col in colnames(scores)) {
        meta_col <- paste0(prefix, col)
        if (meta_col %in% colnames(sobj@meta.data)) {
            log$warn(paste0(
                "Skipping score column '", meta_col,
                "': a metadata column of that name already exists"
            ))
            next
        }
        idx <- match(rownames(scores), cells)
        if (!is.null(subset_cells)) {
            idx[!rownames(scores) %in% subset_cells] <- NA_integer_
        }
        values <- rep(NA_real_, length(cells))
        values[idx[!is.na(idx)]] <- scores[[col]][!is.na(idx)]
        values[!is.finite(values)] <- NA_real_
        sobj@meta.data[[meta_col]] <- values
    }
    sobj
}
### --- scores: end ---

# Apply results to Seurat object in order
case_anno_cols <- list()
case_added_cols <- list()
set_ident_cases <- list()
for (res in results) {
    if (is.null(res)) next

    # Non-DEFAULT cases get their name as the prefix of every column the
    # annotation adds (unless add_prefix is disabled)
    prefix <- ""
    if (!identical(res$name, "DEFAULT") &&
        (is.null(add_prefix) || isTRUE(add_prefix))) {
        prefix <- paste0(res$name, "_")
    }

    cols_before <- colnames(sobj@meta.data)
    sobj <- EmbedSeuratCellTypeAnnotation(
        sobj,
        record = res$record,
        # NULL means "no prefix"; the case name alone carries no other meaning
        case = if (nzchar(prefix)) res$name else NULL,
        ident = res$ident,
        anno_col = res$anno_col,
        # The identity is set at the end instead, where the last enabled case
        # wins
        set_ident = FALSE,
        merge = res$merge,
        subset = res$subset_cells,
        log = log
    )
    # The tool's per-cell scores ride along with the label columns (same
    # prefix, same subset); keep them in the added columns so the report can
    # find the score panel's column
    sobj <- add_case_scores(
        sobj, res$scores, prefix, res$subset_cells, log
    )
    case_added_cols[[res$name]] <- setdiff(
        colnames(sobj@meta.data), cols_before
    )
    case_anno_cols[[res$name]] <- paste0(prefix, res$anno_col)

    # Embed the projections the case's tool produced (mapquery: the reference
    # projection, e.g. ref.umap); cells outside the case's subset keep NA
    for (proj in res$reductions) {
        # With several mapquery cases every reduction a case adds is scoped to
        # that case (`<case>.<reduction>`, the SeuratSubClustering convention),
        # so two cases cannot collide; a name the case already qualified itself
        # (e.g. `CaseA.ref.umap`) is left alone
        name <- proj$name
        if (length(mapquery_cases) > 1 &&
                !startsWith(name, paste0(res$name, "."))) {
            name <- paste0(res$name, ".", name)
        }
        if (name %in% Reductions(sobj)) {
            log$warn(paste0(
                "Skipping projection '", name, "': a reduction of that ",
                "name already exists in the object"
            ))
            next
        }
        emb <- matrix(NA_real_, nrow = ncol(sobj), ncol = ncol(proj$embeddings),
                      dimnames = list(colnames(sobj), colnames(proj$embeddings)))
        # cell names are the ROWS of a projection (columns are the dimensions);
        # a case with no `subset` covers every cell of the payload
        cells <- if (is.null(res$subset_cells)) rownames(proj$embeddings) else res$subset_cells
        if (is.null(cells)) {
            log$warn(paste0(
                "Skipping projection '", name, "': it carries no cell names"
            ))
            next
        }
        idx <- match(cells, rownames(proj$embeddings))
        if (anyNA(idx)) {
            log$warn(paste0(
                "Projection '", name, "': ", sum(is.na(idx)), " of ", length(cells),
                " cell(s) of case '", res$name, "' have no coordinates (left NA)"
            ))
            cells <- cells[!is.na(idx)]
            idx <- idx[!is.na(idx)]
        }
        emb[cells, ] <- proj$embeddings[idx, , drop = FALSE]
        # Keys must stay unique across the object, or two projections would
        # share embedding column names; fall back to a key built from the
        # (already case-qualified) reduction name
        key <- proj$key
        if (key %in% vapply(
            Reductions(sobj), function(r) Key(sobj[[r]]), character(1)
        )) {
            key <- paste0(gsub("[^A-Za-z0-9]", "", toupper(name)), "_")
        }
        sobj[[name]] <- CreateDimReducObject(
            embeddings = emb, key = key, assay = proj$assay
        )
        log$info("Added projection '{name}' ({ncol(sobj)} cells, {sum(!is.na(emb[, 1]))} with coordinates)")
    }

    if (isTRUE(res$set_ident)) {
        set_ident_cases[[length(set_ident_cases) + 1]] <- res$name
    }
}

# Save cluster-to-cell-type mappings from all cluster-based cases
# Collect mappings per case
case_mappings <- list()
all_clusters <- character(0)
for (res in results) {
    if (is.null(res)) next
    if (identical(res$type, "cluster") && !is.null(res$mapping)) {
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
        if (is.null(res) || !identical(res$type, "cluster") || is.null(res$mapping)) next
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
    if (is.null(res)) next

    prefix <- ""
    if (!identical(res$name, "DEFAULT") &&
        (is.null(add_prefix) || isTRUE(add_prefix))) {
        prefix <- paste0(res$name, "_")
    }

    if (identical(res$type, "cluster")) {
        # Both-mode: the tool-specific per-cell column is kept
        if (is.null(res$cells)) next
        cell_case_cols[[res$name]] <- paste0(prefix, colnames(res$cells)[1])
    } else {
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

if (!is.null(post_mutaters) && length(post_mutaters) > 0) {
    log$info("Applying post mutaters ...")
    sobj <- MutateSeuratMeta(sobj, post_mutaters, log = log)
}

log$info("Saving Seurat object...")
save_obj(sobj, outfile)

# Merge the report options over the defaults so a partially-specified env
# still works
report <- list_update(
    list(
        reduction = "dim",
        score_col = NULL,
        purity = TRUE,
        devpars = list(width = 1200, height = 720, res = 100),
        more_formats = "pdf"
    ),
    report
)

# The prefix of the metadata columns a case added
case_prefix <- function(case_name) {
    if (identical(case_name, "DEFAULT") || isFALSE(add_prefix)) {
        return("")
    }
    paste0(case_name, "_")
}

# The per-cell annotation column of a case; NULL when the case only labels
# clusters and carries no per-cell labels
case_label_col <- function(res) {
    if (identical(res$type, "cell")) {
        return(paste0(case_prefix(res$name), res$anno_col))
    }
    if (!is.null(res$cells)) {
        return(paste0(case_prefix(res$name), colnames(res$cells)[1]))
    }
    NULL
}

# The numeric score column a case added: `report$score_col` when given, else
# the first added numeric column (preferring one whose name ends with `.score`)
case_score_col <- function(res) {
    if (!is.null(report$score_col)) {
        if (report$score_col %in% colnames(sobj@meta.data)) {
            return(report$score_col)
        }
        log$warn(paste0(
            "Score column '", report$score_col,
            "' is not a metadata column; detecting the score column instead"
        ))
    }
    added <- setdiff(
        case_added_cols[[res$name]] %||% character(0),
        case_anno_cols[[res$name]]
    )
    added <- intersect(added, colnames(sobj@meta.data))
    numeric_cols <- added[vapply(
        sobj@meta.data[added], is.numeric, logical(1)
    )]
    if (length(numeric_cols) == 0) {
        return(NULL)
    }
    score_cols <- numeric_cols[endsWith(numeric_cols, ".score")]
    if (length(score_cols) > 0) score_cols[[1]] else numeric_cols[[1]]
}

# The purity of a per-cell annotation against the clusters of the identity
# column: the modal label and its share per cluster, with the clusters ordered
# by the factor levels of the identity column (or by first appearance)
compute_purity <- function(labels, clusters) {
    keep <- !is.na(labels) & !is.na(clusters)
    labels <- as.character(labels[keep])
    clusters <- clusters[keep]
    if (length(labels) == 0) {
        return(NULL)
    }
    cluster_levels <- if (is.factor(clusters)) {
        levels(droplevels(clusters))
    } else {
        unique(as.character(clusters))
    }
    rows <- lapply(cluster_levels, function(cluster) {
        in_cluster <- as.character(clusters) == cluster
        counts <- sort(table(labels[in_cluster]), decreasing = TRUE)
        data.frame(
            Cluster = cluster,
            NCells = sum(in_cluster),
            MajorCellType = names(counts)[1],
            NMajor = as.integer(counts[[1]]),
            Purity = as.numeric(counts[[1]]) / sum(in_cluster),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

reporter <- get_reporter()
combine_plots <- getFromNamespace("combine_plots", "plotthis")
report_devpars <- list_update(
    list(width = 1200, height = 720, res = 100),
    report$devpars
)

for (res in results) {
    if (is.null(res)) next

    anno_col <- case_anno_cols[[res$name]]
    if (is.null(anno_col) || !anno_col %in% colnames(sobj@meta.data)) {
        log$warn(paste0(
            "Skipping the report of case '", res$name,
            "': annotation column '", anno_col %||% "<none>",
            "' is not in the metadata"
        ))
        next
    }
    if (all(is.na(sobj@meta.data[[anno_col]]))) {
        log$warn(paste0(
            "Skipping the report of case '", res$name,
            "': annotation column '", anno_col, "' is entirely NA"
        ))
        next
    }

    reduction <- report$reduction
    if (is.null(reduction) || reduction %in% c("dim", "auto")) {
        reduction <- tryCatch(
            SeuratObject::DefaultDimReduc(sobj),
            error = function(e) NULL
        )
    }
    has_reduction <- !is.null(reduction) &&
        isTRUE(reduction %in% Reductions(sobj))
    if (!has_reduction) {
        log$warn(paste0(
            "Skipping the plots of case '", res$name,
            "': reduction '", reduction %||% "<none>",
            "' is not in the object"
        ))
    } else {
        score_col <- case_score_col(res)
        plot_prefix <- file.path(
            outdir, paste0(slugify(res$name), ".annot.umap")
        )
        p_annotation <- scplotter::CellDimPlot(
            sobj,
            group_by = anno_col,
            reduction = reduction,
            title = paste0(res$name, ": annotation")
        )
        if (is.null(score_col)) {
            p <- p_annotation
        } else {
            p_score <- scplotter::FeatureStatPlot(
                sobj,
                features = score_col,
                plot_type = "dim",
                reduction = reduction,
                title = paste0(res$name, ": ", score_col)
            )
            p <- combine_plots(
                list(p_annotation, p_score), combine = TRUE, nrow = 1
            )
        }
        save_plot(
            p, plot_prefix, report_devpars,
            formats = c("png", report$more_formats)
        )

        anno_descr <- paste0(
            "Case `", res$name, "`: the `",
            cases[[res$name]]$tool %||% "unknown",
            "` tool annotated `", anno_col, "`"
        )
        if (!is.null(score_col)) {
            anno_descr <- paste0(
                anno_descr, "; the second panel is coloured by the ",
                "score column `", score_col, "`"
            )
        }
        do.call(
            reporter$add,
            c(
                list(
                    list(kind = "descr", content = paste0(anno_descr, ".")),
                    reporter$image(
                        plot_prefix, report$more_formats, FALSE,
                        kind = "image"
                    )
                ),
                list(h1 = res$name, h2 = "Annotation")
            )
        )
    }

    label_col <- case_label_col(res)
    if (!isTRUE(report$purity) ||
            is.null(label_col) ||
            identical(label_col, default_ident) ||
            is.null(default_ident) ||
            !default_ident %in% colnames(sobj@meta.data)) {
        next
    }
    if (!label_col %in% colnames(sobj@meta.data)) {
        log$warn(paste0(
            "Skipping the purity of case '", res$name,
            "': label column '", label_col,
            "' is not in the metadata (identity column '",
            default_ident, "')"
        ))
        next
    }

    purity <- compute_purity(
        sobj@meta.data[[label_col]],
        sobj@meta.data[[default_ident]]
    )
    if (is.null(purity)) next

    purity_file <- paste0(
        outprefix, ".", slugify(res$name), ".purity.tsv"
    )
    write_table(
        purity, purity_file,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    overall_purity <- sum(purity$NMajor) / sum(purity$NCells)
    log$info(paste0(
        "Purity of case '", res$name, "' against '", default_ident,
        "': ", sprintf("%.3f", overall_purity)
    ))
    reporter$add(
        list(
            kind = "descr",
            content = paste0(
                "Overall purity: ", sprintf("%.3f", overall_purity),
                " (", sum(purity$NCells), " cells, identity column `",
                default_ident, "`)."
            )
        ),
        list(
            kind = "table",
            src = purity_file,
            data = list(nrows = 100)
        ),
        h1 = res$name,
        h2 = "Purity"
    )
}

reporter$save(joboutdir)
