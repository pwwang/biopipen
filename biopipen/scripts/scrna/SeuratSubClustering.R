source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/caching.R")

library(Seurat)
library(future)
library(rlang)
library(tidyr)
library(dplyr)
library(tidyseurat)
library(digest)

set.seed(8525)

srtfile <- {{in.srtobj | quote}}
rdsfile <- {{out.rdsfile | quote}}
joboutdir <- {{job.outdir | quote}}
envs <- {{envs | r: todot = "-"}}

options(str = strOptions(vec.len = 5, digits.d = 5))
options(future.globals.maxSize = 80000 * 1024^2)
plan(strategy = "multicore", workers = envs$ncores)

.expand_dims <- function(args, name = "dims") {
    # Expand dims from 30 to 1:30
    if (!is.null(args) && is.numeric(args[[name]]) && length(args[[name]] == 1)) {
        args[[name]] <- 1:args[[name]]
    }
    args
}

envs$RunUMAP <- .expand_dims(envs$RunUMAP)
envs$FindNeighbors <- .expand_dims(envs$FindNeighbors)

log_info("Reading Seurat object ...")
srtobj <- readRDS(srtfile)

if (isTRUE(envs$cache)) { envs$cache <- joboutdir }
if (length(envs$cache) > 1) {
    log_warn("Multiple cache directories (envs.cache) detected, using the first one.")
    envs$cache <- envs$cache[1]
}

if (!is.null(envs$mutaters) && length(envs$mutaters) > 0) {
    log_info("Mutating Seurat object ...")
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

if (length(envs$cases) == 0) {
    envs$cases <- list(subcluster = list())
}

for (key in names(envs$cases)) {
    log_info("")
    log_info("Running case '{key}' ...")
    log_info("===========================================")
    case <- envs$cases[[key]]
    case$RunUMAP <- .expand_dims(case$RunUMAP)
    case$FindNeighbors <- .expand_dims(case$FindNeighbors)

    case <- list_update(
        list(
            subset = envs$subset,
            RunUMAP = envs$RunUMAP,
            FindNeighbors = envs$FindNeighbors,
            FindClusters = envs$FindClusters
        ),
        case
    )

    if (is.null(case$subset) || length(case$subset) == 0) {
        stop(paste0("`subset` for case '", key, "' is empty."))
    }

    log_info("- Subsetting ...")
    sobj <- tryCatch({
        srtobj %>% filter(!!parse_expr(case$subset))
    }, error = function(e) {
        stop(paste0("  Error in subset: ", e$message))
    })
    sobj_sig <- capture.output(str(sobj))
    dig_sig <- digest::digest(sobj_sig, algo = "md5")
    dig_sig <- substr(dig_sig, 1, 8)
    cache_dir <- NULL
    if (is.character(envs$cache)) {
        cache_dir <- file.path(envs$cache, paste0(dig_sig, ".seurat_cache"))
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        writeLines(sobj_sig, file.path(cache_dir, "signature.txt"))
    }

    cached <- get_cached(case$RunUMAP, "RunUMAP", cache_dir)
    reduc_name <- case$RunUMAP$reduction.name %||% "umap"
    if (is.null(cached$data)) {
        log_info("- Running RunUMAP ...")
        umap_args <- list_setdefault(
            case$RunUMAP,
            object = sobj,
            dims = 1:30,
            reduction = sobj@misc$integrated_new_reduction %||% "pca"
        )
        ncells <- ncol(sobj)
        umap_args$dims <- 1:min(max(umap_args$dims), ncells - 1)
        umap_method <- case$RunUMAP$umap.method %||% "uwot"
        if (umap_method == "uwot" && is.null(case$RunUMAP$n.neighbors)) {
            # https://github.com/satijalab/seurat/issues/4312
            umap_args$n.neighbors <- min(ncells - 1, 30)
        }
        sobj <- do_call(RunUMAP, umap_args)
        cached$data <- list(reduc = sobj@reductions[[reduc_name]], commands = sobj@commands)
        save_to_cache(cached, "RunUMAP", cache_dir)
    } else {
        log_info("- Loading cached RunUMAP ...")
        sobj@reductions[[reduc_name]] <- cached$data$reduc
        sobj@commands <- cached$data$commands
    }
    reduc <- cached$data$reduc

    cached <- get_cached(case$FindNeighbors, "FindNeighbors", cache_dir)
    if (is.null(cached$data)) {
        log_info("- Running FindNeighbors ...")
        case$FindNeighbors$object <- sobj
        if (is.null(case$FindNeighbors$reduction)) {
            case$FindNeighbors$reduction <- sobj@misc$integrated_new_reduction %||% "pca"
        }
        sobj <- do_call(FindNeighbors, case$FindNeighbors)
        cached$data <- list(graphs = sobj@graphs, commands = sobj@commands)
        save_to_cache(cached, "FindNeighbors", cache_dir)
    } else {
        log_info("- Loading cached FindNeighbors ...")
        sobj@graphs <- cached$data$graphs
        sobj@commands <- cached$data$commands
    }

    case$FindClusters$random.seed <- case$FindClusters$random.seed %||% 8525
    resolution <- case$FindClusters$resolution %||% 0.8
    if (is.character(resolution)) {
        if (grepl(",", resolution)) {
            resolution <- as.numeric(trimws(unlist(strsplit(resolution, ","))))
        } else {
            resolution <- as.numeric(resolution)
        }
    }
    for (res in resolution) {
        case$FindClusters$resolution <- res
        cached <- get_cached(case$FindClusters, paste0("FindClusters_", res), cache_dir)
        res_key <- paste0("seurat_clusters_", res)
        if (is.null(cached$data)) {
            log_info("- Running FindClusters at resolution: {res} ...")
            case$FindClusters$object <- sobj
            sobj1 <- do_call(FindClusters, case$FindClusters)
            levels(sobj1$seurat_clusters) <- paste0("s", as.numeric(levels(sobj1$seurat_clusters)) + 1)
            sobj1[[res_key]] <- sobj1$seurat_clusters
            cached$data <- sobj1@meta.data[, res_key, drop = FALSE]
            save_to_cache(cached, paste0("FindClusters_", res), cache_dir)
        } else {
            log_info("- Using cached FindClusters at resolution: {res} ...")
        }
        ident_table <- table(cached$data[[res_key]])
        log_info("  Found {length(ident_table)} clusters")
        print(ident_table)
        cat("\n")
    }
    log_info("- Updating meta.data with subclusters...")
    srtobj <- AddMetaData(srtobj, metadata = cached$data, col.name = key)
    srtobj[[paste0("sub_umap_", key)]] <- reduc
}

log_info("Saving results ...")
saveRDS(srtobj, file = rdsfile)
