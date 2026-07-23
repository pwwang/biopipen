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
backup_col <- {{envs.backup_col | r}}
merge <- {{envs.merge | r}}
newcol <- {{envs.newcol | r}}
add_prefix <- {{envs.add_prefix | r}}
outtype <- {{envs.outtype | r}}
cases <- {{envs.cases | r}}
ncores <- {{envs.ncores | r}}

assay <- {{envs.assay | r}}
# Tool-specific defaults (Jinja2-substituted)
sctype_tissue <- {{envs.sctype_tissue | r}}
sctype_db <- {{envs.sctype_db | r}}
hitype_tissue <- {{envs.hitype_tissue | r}}
hitype_db <- {{envs.hitype_db | r}}
scsorter_db <- {{envs.scsorter_db | r}}
scsorter_args <- {{envs.scsorter_args | r}}
sccatch_args <- {{envs.sccatch_args | r}}
celltypist_args <- {{envs.celltypist_args | r}}
scina_db <- {{envs.scina_db | r}}
scina_args <- {{envs.scina_args | r}}
cell_types <- {{envs.cell_types | r}}
more_cell_types <- {{envs.more_cell_types | r}}

log <- get_logger()

# Source all tool function definitions
biopipen_dir <- {{ biopipen_dir | r }}
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-hitype.R"))
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-sctype.R"))
source(file.path(biopipen_dir, "scripts", "scrna", "sctype.R"))
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-sccatch.R"))
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-celltypist.R"))
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-scsorter.R"))
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-direct.R"))
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-cell.R"))
source(file.path(biopipen_dir, "scripts", "scrna", "CellTypeAnnotation-scina.R"))

# Build defaults list for expand_cases (exclude ncores — not inherited)
defaults <- list(
    tool = tool,
    ident = ident,
    backup_col = backup_col,
    merge = merge,
    newcol = newcol,
    sctype_tissue = sctype_tissue,
    sctype_db = sctype_db,
    hitype_tissue = hitype_tissue,
    hitype_db = hitype_db,
    scsorter_db = scsorter_db,
    scsorter_args = scsorter_args,
    sccatch_args = sccatch_args,
    scina_db = scina_db,
    scina_args = scina_args,
    celltypist_args = celltypist_args,
    cell_types = cell_types,
    more_cell_types = more_cell_types
)

cases <- expand_cases(cases, defaults, default_case = "DEFAULT")

# Cluster-based tools
CLUSTER_TOOLS <- c("hitype", "sctype", "sccatch", "scina", "direct")

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

outdir <- dirname(outfile)
outprefix <- file.path(outdir, tools::file_path_sans_ext(basename(outfile)))

# Pre-convert to h5ad if any case needs it
needs_h5ad <- any(sapply(cases, function(c)
    identical(c$tool, "celltypist")
))
h5ad_path <- NULL
if (needs_h5ad) {
    log$info("Pre-converting Seurat object to h5ad for Python-based tools...")
    h5ad_path <- paste0(outprefix, ".shared.h5ad")
    ConvertSeuratToAnnData(
        sobj,
        outfile = h5ad_path,
        assay = celltypist_args$assay,
        log = log
    )
}

# Resolve ident for each case and determine tool type
cases <- lapply(cases, function(case) {
    if (is.null(case$ident)) {
        case$ident <- GetIdentityColumn(sobj)
    }
    case$is_cluster_based <- case$tool %in% CLUSTER_TOOLS ||
        (identical(case$tool, "celltypist") &&
         (!is.null(case$celltypist_args$over_clustering) &&
          !isFALSE(case$celltypist_args$over_clustering)))
    case
})

# Run each case
run_case <- function(case_name) {
    case <- cases[[case_name]]

    log$info("Running annotation for case: {case_name}")
    case$celltypist_args$assay <- case$celltypist_args$assay %||% assay
    case$scsorter_args$assay <- case$scsorter_args$assay %||% assay

    result <- switch(case$tool,
        hitype = annotate_hitype(
            sobj, case$ident, case$hitype_tissue, case$hitype_db
        ),
        sctype = annotate_sctype(
            sobj, case$ident, case$sctype_tissue, case$sctype_db
        ),
        sccatch = annotate_sccatch(
            sobj, case$ident, case$sccatch_args
        ),
        celltypist = annotate_celltypist(
            sobj, case$ident, case$celltypist_args, outdir,
            h5ad_path = h5ad_path, case_id = case_name
        ),
        scsorter = annotate_scsorter(
            sobj, case$ident, case$scsorter_db, case$scsorter_args
        ),
        scina = annotate_scina(
            sobj, case$ident, case$scina_db, case$scina_args
        ),
        direct = annotate_direct(
            sobj, case$ident, case$cell_types, case$more_cell_types
        ),
        cell = annotate_cell(sobj, case$cell_types),
        stop(paste0("Unknown tool: ", case$tool))
    )

    list(
        name = case_name,
        mapping = result$mapping,
        cell_annotations = result$cell_annotations,
        annotation_col = result$annotation_col,
        more = result$more,
        ident = case$ident,
        newcol = case$newcol,
        merge = case$merge,
        backup_col = case$backup_col,
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
last_annot_col <- NULL
last_annot_ident <- NULL
for (res in results) {
    if (is.null(res)) next

    # Determine prefix for non-DEFAULT cases
    prefix <- ""
    if (!identical(res$name, "DEFAULT") &&
        (is.null(add_prefix) || isTRUE(add_prefix))) {
        prefix <- paste0(res$name, "_")
    }

    if (res$is_cluster_based) {
        # Cluster-based tools
        save_as_column <- NULL
        if (!is.null(res$newcol)) {
            save_as_column <- paste0(prefix, res$newcol)
            log$info(
                "Adding annotation as new column: {save_as_column}"
            )
            # Use backup_col with prefix for non-DEFAULT, but only for
            # the first case to avoid overwriting the backup
            bcol <- res$backup_col
            if (!is.null(bcol)) {
                bcol <- paste0(prefix, bcol)
            }
            sobj <- RenameSeuratIdents(
                sobj,
                mapping = res$mapping,
                ident = res$ident,
                save_as = save_as_column,
                merge = res$merge,
                backup = bcol
            )
        } else {
            # No newcol: overwrite ident values in-place
            log$info(
                "Overwriting ident column '{res$ident}' with annotations"
            )
            bcol <- res$backup_col
            if (!is.null(bcol)) {
                bcol <- paste0(prefix, bcol)
            }
            sobj <- RenameSeuratIdents(
                sobj,
                mapping = res$mapping,
                ident = res$ident,
                merge = res$merge,
                backup = bcol
            )
            save_as_column <- res$ident
        }
        last_annot_col <- save_as_column
        last_annot_ident <- res$ident
    } else {
        # Cell-based tools
        if (!is.null(res$cell_annotations)) {
            cell_annotations <- res$cell_annotations
            # Prefix column names for non-DEFAULT cases
            if (nzchar(prefix)) {
                colnames(cell_annotations) <- paste0(
                    prefix, colnames(cell_annotations)
                )
                annotation_col <- paste0(prefix, res$annotation_col)
            } else {
                annotation_col <- res$annotation_col
            }

            # Add to meta.data
            log$info(
                "Adding cell-level annotations to metadata"
            )
            # Remove existing columns that will be overwritten
            existing <- intersect(
                colnames(cell_annotations),
                colnames(sobj@meta.data)
            )
            if (length(existing) > 0) {
                sobj@meta.data[, existing] <- NULL
            }
            sobj@meta.data <- cbind(sobj@meta.data, cell_annotations)

            if (!is.null(res$newcol)) {
                # Copy annotation to newcol
                save_as_column <- paste0(prefix, res$newcol)
                sobj@meta.data[[save_as_column]] <- sobj@meta.data[[annotation_col]]
                log$info(
                    "Copied annotation to new column: {save_as_column}"
                )
            } else {
                save_as_column <- annotation_col
            }
            last_annot_col <- save_as_column
            last_annot_ident <- NULL
        }
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
        }
    }
}

# Set identity to last case's annotation column
if (!is.null(last_annot_col)) {
    log$info("Setting identity to last case annotation: {last_annot_col}")
    Idents(sobj) <- last_annot_col
}

log$info("Saving Seurat object...")
save_obj(sobj, outfile)
