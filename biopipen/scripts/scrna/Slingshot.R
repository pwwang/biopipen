library(rlang)
library(Seurat)
library(slingshot)
library(biopipen.utils)

sobjfile <- {{in.sobjfile | r}}
outfile <- {{out.outfile | r}}
group_by <- {{envs.group_by | r}}
reduction <- {{envs.reduction | r}}
dims <- {{envs.dims | r}}
start <- {{envs.start | r}}
end <- {{envs.end | r}}
cases <- {{envs.cases | r}}
reverse <- {{envs.reverse | r}}
align_start <- {{envs.align_start | r}}
subset <- {{envs.subset | r}}
split_by <- {{envs.split_by | r}}
seed <- {{envs.seed | r}}
ncores <- {{envs.ncores | r}}

set.seed(seed)
qs2::qopt("nthreads", value = ncores)

log <- get_logger()

log$info("Reading Seurat object ...")
srtobj <- read_obj(sobjfile)
ident <- GetIdentityColumn(srtobj)

defaults <- list(
    group_by = group_by,
    reduction = reduction,
    dims = dims,
    start = start,
    end = end,
    reverse = reverse,
    align_start = align_start,
    subset = subset,
    split_by = split_by,
    seed = seed
)

cases <- expand_cases(cases, defaults, default_case = "DEFAULT")


do_case_split <- function(prefix, case, srt_split, split = NULL) {
    clusterLabels <- as.character(srt_split[[case$group_by, drop = TRUE]])
    small_clusters <- names(table(clusterLabels))[table(clusterLabels) < 2]
    if (length(small_clusters) > 0) {
        log$warn("- Dropping clusters with less than 2 cells: {paste(small_clusters, collapse = ', ')}")
        srt_split <- scplotter:::subset_seurat(srt_split, subset = !(!!rlang::sym(case$group_by) %in% small_clusters))
    }

    if (!is.null(split)) {
        log$info("  Running Slingshot for split: {split}")
    } else {
        log$info("  Running Slingshot ...")
    }

    sl <- slingshot(
        data = if (!is.null(case$dims)) {
            as.data.frame(srt_split[[case$reduction]]@cell.embeddings[, case$dims])
        } else {
            as.data.frame(srt_split[[case$reduction]]@cell.embeddings)
        },
        clusterLabels = as.character(srt_split[[case$group_by, drop = TRUE]]),
        start.clus = case$start, end.clus = case$end
    )

    df <- as.data.frame(slingPseudotime(sl))
    colnames(df) <- paste0(prefix, colnames(df))
    log$info("  Lineages inferred: {paste(colnames(df), collapse = ', ')}")
    if (isTRUE(case$reverse)) {
        if (isTRUE(case$align_start)) {
            df <- apply(df, 2, function(x) max(x, na.rm = TRUE) - x)
        } else {
            df <- max(df, na.rm = TRUE) - df
        }
    }

    df[[paste0(prefix, "BranchID")]] <- slingBranchID(sl)
    return(df)
}

do_case <- function(prefix) {
    log$info("- Running case: {prefix}")
    case <- cases[[prefix]]
    if (identical(prefix, "DEFAULT")) {
        prefix <- ""
    }
    if (nzchar(prefix)) {
        prefix <- paste0(prefix, "_")
    }

    case$group_by <- case$group_by %||% ident
    if (!case$group_by %in% colnames(srtobj@meta.data)) {
        stop(paste("Grouping column", case$group_by, "not found in the Seurat object"))
    } else {
        log$info("  Using group_by: {case$group_by}")
    }
    case$reduction <- case$reduction %||% scplotter:::default_dimreduc(srtobj)
    if (!is.null(case$dims)) {
        case$dims <- biopipen.utils:::.expand_number(case$dims)
    }

    log$info("  Filtering cells in NA group_by ...")
    srt_sub <- srtobj[, !is.na(srtobj[[case$group_by, drop = TRUE]])]
    if (!is.null(case$subset)) {
        srt_sub <- scplotter:::subset_seurat(srt_sub, subset = !!rlang::parse_expr(case$subset))
    }

    if (!is.null(case$split_by)) {
        result <- NULL
        if (!case$split_by %in% colnames(srt_sub@meta.data)) {
            stop(paste("Split column", case$split_by, "not found in the Seurat object"))
        }
        splits <- unique(srt_sub@meta.data[[case$split_by]])
        splits <- splits[!is.na(splits)]
        for (split in splits) {
            srt_split <- scplotter:::subset_seurat(srt_sub, subset = !!rlang::sym(case$split_by) == split)
            result_split <- do_case_split(prefix, case, srt_split, split = split)
            # result <- rbind(result, result_split)
            # <prefix>_Lineage1, <prefix>_Lineage2, ..., <prefix>_BranchID
            # The number of lineages may vary for different splits, so we need to handle this case.
            if (is.null(result)) {
                result <- result_split
            } else {
                # Align columns by name, fill missing columns with NA
                all_cols <- union(colnames(result), colnames(result_split))
                missing_cols_result <- setdiff(all_cols, colnames(result))
                missing_cols_split <- setdiff(all_cols, colnames(result_split))
                if (length(missing_cols_result) > 0) {
                    result[missing_cols_result] <- NA
                }
                if (length(missing_cols_split) > 0) {
                    result_split[missing_cols_split] <- NA
                }
                result <- rbind(result, result_split)
            }
        }
        return(result)
    } else {
        do_case_split(prefix, case, srt_sub)
    }
}

for (prefix in names(cases)) {
    result <- do_case(prefix)

    log$info("  Attaching result to Seurat object ...")
    srtobj <- AddMetaData(srtobj, metadata = result)
}

srtobj <- AddSeuratCommand(srtobj, "Slingshot", "slingshot(...)")

log$info("Saving Seurat object ...")
save_obj(srtobj, outfile)
