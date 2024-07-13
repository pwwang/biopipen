{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "caching.R" | source_r }}

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

.expand_resolution <- function(resolution) {
    expanded_res <- c()
    for (res in resolution) {
        if (is.numeric(res)) {
            expanded_res <- c(expanded_res, res)
        } else {
            # is.character
            parts <- trimws(unlist(strsplit(res, ",")))
            for (part in parts) {
                if (grepl(":", part)) {
                    ps <- trimws(unlist(strsplit(part, ":")))
                    if (length(ps) == 2) { ps <- c(ps, 0.1) }
                    if (length(ps) != 3) {
                        stop("Invalid resolution format: {part}. Expected 2 or 3 parts separated by ':' for a range.")
                    }
                    ps <- as.numeric(ps)
                    expanded_res <- c(expanded_res, seq(ps[1], ps[2], by = ps[3]))
                } else {
                    expanded_res <- c(expanded_res, as.numeric(part))
                }
            }
        }
    }
    # keep the last resolution at last
    rev(unique(rev(round(expanded_res, 2))))
}

# recode clusters from 0, 1, 2, ... to s1, s2, s3, ...
.recode_clusters <- function(clusters) {
    recode <- function(x) paste0("s", as.integer(as.character(x)) + 1)
    clusters <- factor(recode(clusters), levels = recode(levels(clusters)))
    clusters
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
    resolution <- case$FindClusters$resolution <- .expand_resolution(case$FindClusters$resolution %||% 0.8)
    cached <- get_cached(case$FindClusters, "SubClustering", cache_dir)
    if (is.null(cached$data)) {
        log_info("- Running FindClusters at resolution: {paste(resolution, collapse = ',')} ...")
        case$FindClusters$object <- sobj
        case$FindClusters$cluster.name <- paste0(key, ".", resolution)
        # use sobj1 to avoid overwriting the previous clustering results (as they have the same graph name
        sobj1 <- do_call(FindClusters, case$FindClusters)
        sobj1@meta.data[[key]] <- .recode_clusters(sobj1@meta.data$seurat_clusters)
        for (clname in case$FindClusters$cluster.name) {
            sobj1@meta.data[[clname]] <- .recode_clusters(sobj1@meta.data[[clname]])
        }
        cached$data <- list(
            clusters = sobj1@meta.data[, c(case$FindClusters$cluster.name, key), drop = FALSE],
            command = sobj1@commands$FindClusters
        )
        save_to_cache(cached, "SubClustering", cache_dir)
        rm(sobj1)
    } else {
        log_info("- Using cached FindClusters at resolution: {paste(resolution, collapse = ',')} ...")
    }

    ident_table <- table(cached$data$clusters[[key]])
    log_info("  Found {length(ident_table)} clusters at resolution: {resolution[length(resolution)]}")
    print(ident_table)
    cat("\n")

    log_info("- Updating meta.data with subclusters...")
    srtobj <- AddMetaData(srtobj, metadata = cached$data$clusters)
    srtobj[[paste0("sub_umap_", key)]] <- reduc
    srtobj@commands[[paste0("FindClusters.", key)]] <- cached$data$command
}

log_info("Saving results ...")
saveRDS(srtobj, file = rdsfile)
