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
seed <- {{envs.seed | r}}

set.seed(seed)

log <- get_logger()

log$info("Reading Seurat object ...")
srt <- read_obj(sobjfile)

defaults <- list(
    group_by = group_by,
    reduction = reduction,
    dims = dims,
    start = start,
    end = end,
    reverse = reverse,
    align_start = align_start,
    seed = seed
)

cases <- expand_cases(cases, defaults, default_case = "DEFAULT")

do_case <- function(prefix) {
    log$info("Running Slingshot for case: {prefix}")
    case <- cases[[prefix]]
    if (identical(prefix, "DEFAULT")) {
        prefix <- ""
    }
    if (nzchar(prefix)) {
        prefix <- paste0(prefix, "_")
    }

    case$group_by <- case$group_by %||% GetIdentityColumn(srt)
    if (!case$group_by %in% colnames(srt@meta.data)) {
        stop(paste("Grouping column", case$group_by, "not found in the Seurat object"))
    } else {
        log$info("- Using group_by: {case$group_by}")
    }
    case$reduction <- case$reduction %||% scplotter:::default_dimreduc(srt)
    if (!is.null(case$dims)) {
        case$dims <- biopipen.utils:::.expand_number(case$dims)
    }

    log$info("- Filtering cells in NA group_by ...")
    srt_sub <- srt[, !is.na(srt[[case$group_by, drop = TRUE]])]

    clusterLabels <- as.character(srt_sub[[case$group_by, drop = TRUE]])
    small_clusters <- names(table(clusterLabels))[table(clusterLabels) < 2]
    if (length(small_clusters) > 0) {
        log$warn("- Dropping clusters with less than 2 cells: {paste(small_clusters, collapse = ', ')}")
        srt_sub <- scplotter:::subset_seurat(srt_sub, subset = !(!!rlang::sym(case$group_by) %in% small_clusters))
    }

    log$info("- Running Slingshot ...")
    sl <- slingshot(
        data = if (!is.null(case$dims)) {
            as.data.frame(srt_sub[[case$reduction]]@cell.embeddings[, case$dims])
        } else {
            as.data.frame(srt_sub[[case$reduction]]@cell.embeddings)
        },
        clusterLabels = as.character(srt_sub[[case$group_by, drop = TRUE]]),
        start.clus = case$start, end.clus = case$end
    )

    df <- as.data.frame(slingPseudotime(sl))
    colnames(df) <- paste0(prefix, colnames(df))
    log$info("- Lineages inferred: {paste(colnames(df), collapse = ', ')}")
    if (isTRUE(case$reverse)) {
        if (isTRUE(case$align_start)) {
            df <- apply(df, 2, function(x) max(x, na.rm = TRUE) - x)
        } else {
            df <- max(df, na.rm = TRUE) - df
        }
    }

    list(df = df, branch = slingBranchID(sl), branch_name = paste0(prefix, "BranchID"))
}

for (prefix in names(cases)) {
    result <- do_case(prefix)

    log$info("- Attaching result to Seurat object ...")
    srt <- AddMetaData(srt, metadata = result$df)
    srt <- AddMetaData(srt, metadata = result$branch, col.name = result$branch_name)
}

srt <- AddSeuratCommand(srt, "Slingshot", "slingshot(...)")

log$info("Saving Seurat object ...")
save_obj(srt, outfile)
