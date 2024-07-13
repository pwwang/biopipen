
expand_dims <- function(args, name = "dims") {
    # Expand dims from 30 to 1:30
    if (is.numeric(args[[name]]) && length(args[[name]] == 1)) {
        args[[name]] <- 1:args[[name]]
    }
    args
}

expand_resolution <- function(resolution) {
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

# recode clusters from 0, 1, 2, ... to c1, c2, c3, ...
recode_clusters <- function(clusters) {
    recode <- function(x) paste0("c", as.integer(as.character(x)) + 1)
    clusters <- factor(recode(clusters), levels = recode(levels(clusters)))
    clusters
}

run_transformation <- function(sobj) {
    if (length(envs$ScaleData) == 0 && length(envs$SCTransform) == 0) {
        log_warn("Skipping ScaleData/SCTransform (neither specified) ...")
        return(sobj)
    }
    if (length(envs$ScaleData) > 0 && length(envs$SCTransform) > 0) {
        stop("Both envs.ScaleData and envs.SCTransform are specified. Please choose either.")
    }
    if (length(envs$ScaleData) > 0) {
        if (DefaultAssay(sobj) == "SCT") {
            stop("SCT assay detected, but envs.ScaleData is specified. Use envs.SCTransform instead.")
        }
        cached <- get_cached(envs$ScaleData, "ScaleData", cache_dir)
        if (is.null(cached$data)) {
            log_info("Running ScaleData ...")
            sobj <- do_call(ScaleData, c(list(object = sobj), envs$ScaleData))
            cached$data <- list(assay = sobj@assays$RNA, commands = sobj@commands)
            save_to_cache(cached, "ScaleData", cache_dir)
        } else {
            log_info("Loading cached ScaleData ...")
            sobj@assays$RNA <- cached$data$assay
            sobj@commands <- cached$data$commands
            DefaultAssay(sobj) <- "RNA"
        }
    } else if (length(envs$SCTransform) > 0) {
        if (DefaultAssay(sobj) != "SCT") {
            stop("SCT assay not detected, but envs.SCTransform is specified. Use envs.ScaleData instead.")
        }
        cached <- get_cached(envs$SCTransform, "SCTransform", cache_dir)
        asssay <- envs$SCTransform$new.assay.name %||% "SCT"
        if (is.null(cached$data)) {
            log_info("Running SCTransform ...")
            sobj <- do_call(SCTransform, c(list(object = sobj), envs$SCTransform))
            cached$data <- list(assay = sobj@assays$SCT, commands = sobj@commands)
            save_to_cache(cached, "SCTransform", cache_dir)
        } else {
            log_info("Loading cached SCTransform ...")
            sobj@assays[[assay]] <- cached$data$assay
            sobj@commands <- cached$data$commands
            DefaultAssay(sobj) <- assay
        }
    }
    sobj
}

run_umap <- function(sobj) {
    cached <- get_cached(
        list(sobj = sobj, RunUMAP = envs$RunUMAP),
        "RunUMAP",
        cache_dir
    )
    reduc_name <- envs$RunUMAP$reduction.name %||% "umap"
    if (is.null(cached$data)) {
        log_info("Running RunUMAP ...")
        umap_args <- list_setdefault(
            envs$RunUMAP,
            object = sobj,
            dims = 1:30,
            reduction = sobj@misc$integrated_new_reduction %||% "pca"
        )
        ncells <- ncol(sobj)
        umap_args$dims <- 1:min(max(umap_args$dims), ncells - 1)
        umap_method <- envs$RunUMAP$umap.method %||% "uwot"
        if (umap_method == "uwot" && is.null(envs$RunUMAP$n.neighbors)) {
            # https://github.com/satijalab/seurat/issues/4312
            umap_args$n.neighbors <- min(ncells - 1, 30)
        }
        sobj <- do_call(RunUMAP, umap_args)
        cached$data <- list(reduc = sobj@reductions[[reduc_name]], commands = sobj@commands)
        save_to_cache(cached, "RunUMAP", cache_dir)
    } else {
        log_info("Loading cached RunUMAP ...")
        sobj@reductions[[reduc_name]] <- cached$data$reduc
        sobj@commands <- cached$data$commands
    }

    sobj
}

run_findneighbors <- function(sobj) {
    cached <- get_cached(
        list(sobj = sobj, FindNeighbors = envs$FindNeighbors),
        "FindNeighbors",
        cache_dir
    )
    if (is.null(cached$data)) {
        log_info("Running FindNeighbors ...")
        envs$FindNeighbors$object <- sobj
        envs$FindNeighbors$reduction <- sobj@misc$integrated_new_reduction %||% "pca"
        sobj <- do_call(FindNeighbors, envs$FindNeighbors)
        cached$data <- list(graphs = sobj@graphs, commands = sobj@commands)
        save_to_cache(cached, "FindNeighbors", cache_dir)
    } else {
        log_info("Loading cached FindNeighbors ...")
        sobj@graphs <- cached$data$graphs
        sobj@commands <- cached$data$commands
    }

    sobj
}

run_findclusters <- function(sobj) {
    cached <- get_cached(
        list(sobj = sobj, FindClusters = envs$FindClusters),
        "FindClusters",
        cache_dir
    )
    if (is.null(cached$data)) {
        findclusters_args <- envs$FindClusters
        findclusters_args$random.seed <- findclusters_args$random.seed %||% 8525
        resolution <- findclusters_args$resolution <- expand_resolution(findclusters_args$resolution %||% 0.8)
        log_info("Running FindClusters at resolution: {paste(resolution, collapse=',')} ...")

        findclusters_args$object <- sobj
        findclusters_args$cluster.name <- paste0("seurat_clusters.", resolution)
        sobj <- do_call(FindClusters, findclusters_args)

        for (clname in findclusters_args$cluster.name) {
            sobj@meta.data[[clname]] <- recode_clusters(sobj@meta.data[[clname]])
        }
        sobj@meta.data$seurat_clusters <- recode_clusters(sobj@meta.data$seurat_clusters)
        Idents(sobj) <- "seurat_clusters"

        ident_table <- table(Idents(sobj))
        log_info("- Found {length(ident_table)} clusters at resolution {resolution[length(resolution)]}")
        print(ident_table)
        cat("\n")

        cached$data <- list(
            clusters = sobj@meta.data[, c(findclusters_args$cluster.name, "seurat_clusters"), drop = FALSE],
            commands = sobj@commands
        )
        save_to_cache(cached, "FindClusters", cache_dir)
    } else {
        log_info("Loading cached FindClusters ...")

        sobj <- AddMetaData(sobj, metadata = cached$data$clusters)
        Idents(sobj) <- "seurat_clusters"
        sobj@commands <- cached$data$commands
    }

    sobj
}

run_prepsctfindmarkers <- function(sobj) {
    if (DefaultAssay(sobj) == "SCT") {
        cached <- get_cached(list(sobj = sobj), "PrepSCTFindMarkers", cache_dir)
        if (is.null(cached$data)) {
            # https://github.com/satijalab/seurat/issues/6968
            log_info("Running PrepSCTFindMarkers ...")
            sobj <- PrepSCTFindMarkers(sobj)
            # compose a new SeuratCommand to record it to sobj@commands
            scommand <- sobj@commands$FindClusters
            scommand@name <- "PrepSCTFindMarkers"
            scommand@time.stamp <- Sys.time()
            scommand@assay.used <- "SCT"
            scommand@call.string <- "PrepSCTFindMarkers(object = sobj)"
            scommand@params <- list()
            sobj@commands$PrepSCTFindMarkers <- scommand

            cached$data <- sobj
            save_to_cache(cached, "PrepSCTFindMarkers", cache_dir)
        } else {
            log_info("Loading cached PrepSCTFindMarkers ...")
            sobj <- cached$data
        }
    }

    sobj
}
