
stringify_list <- function(x) {
    paste(sapply(names(x), function(n) paste(n, x[[n]], sep = " = ") ), collapse = "; ")
}

format_args <- function(args) {
    paste(capture.output(str(args)), collapse = ", ")
}

rename_files = function(e, sample, path) {
    tmpdatadir = file.path(joboutdir, "renamed", sample)
    if (dir.exists(tmpdatadir)) {
        unlink(tmpdatadir, recursive = TRUE)
    }
    dir.create(tmpdatadir, recursive = TRUE, showWarnings = FALSE)
    barcodefile = Sys.glob(file.path(path, "*barcodes.tsv.gz"))[1]
    file.symlink(
        normalizePath(barcodefile),
        file.path(tmpdatadir, "barcodes.tsv.gz")
    )
    genefile = glob(file.path(path, "*{genes,features}.tsv.gz"))[1]
    file.symlink(
        normalizePath(genefile),
        file.path(tmpdatadir, "features.tsv.gz")
    )
    matrixfile = Sys.glob(file.path(path, "*matrix.mtx.gz"))[1]
    file.symlink(
        normalizePath(matrixfile),
        file.path(tmpdatadir, "matrix.mtx.gz")
    )
    Read10X(data.dir = tmpdatadir)
}


perform_cell_qc <- function(sobj, per_sample = FALSE) {
    log_prefix <- ifelse(per_sample, "  ", "- ")
    log_info("{log_prefix}Adding metadata for QC ...")
    sobj$percent.mt <- PercentageFeatureSet(sobj, pattern = "^MT-|^Mt-|^mt-")
    sobj$percent.ribo <- PercentageFeatureSet(sobj, pattern = "^RP[SL]|^Rp[sl]")
    sobj$percent.hb <- PercentageFeatureSet(sobj, pattern = "^HB[^P]|^Hb[^p]")
    sobj$percent.plat <- PercentageFeatureSet(sobj, pattern = "PECAM1|PF4|Pecam1|Pf4")

    if (is.null(envs$cell_qc) || length(envs$cell_qc) == 0) {
        log_warn("{log_prefix}No cell QC criteria is provided. All cells will be kept.")
        cell_qc <- "TRUE"
    } else {
        cell_qc <- envs$cell_qc
    }

    sobj@meta.data <- sobj@meta.data %>% mutate(.QC = !!rlang::parse_expr(cell_qc))

    if (is.null(cell_qc_df)) {
        cell_qc_df <<- sobj@meta.data[, c("Sample", ".QC", feats), drop = FALSE]
    } else {
        cell_qc_df <<- rbind(cell_qc_df, sobj@meta.data[, c("Sample", ".QC", feats), drop = FALSE])
    }

    # Do the filtering
    log_info("{log_prefix}Filtering cells using QC criteria ...")
    sobj <- subset(sobj, subset = .QC)
    sobj$.QC <- NULL

    return(sobj)
}

report_cell_qc = function(ngenes) {
    # uses cell_qc_df

    # Violin plots
    log_info("- Plotting violin plots ...")
    add_report(
        list(
            kind = "descr",
            content = paste(
                "The violin plots for each feature. The cells are grouped by sample.",
                "The cells that fail the QC criteria are colored in red, and",
                "the cells that pass the QC criteria are colored in black.",
                "The cells that fail the QC criteria are filtered out in the returned Seurat object."
            )
        ),
        h1 = "Violin Plots"
    )
    for (feat in feats) {
        log_info("  For feature: {feat}")
        vln_p <- ggplot(cell_qc_df, aes(x = Sample, y = !!sym(feat), color = .QC)) +
            geom_violin(fill = "white", width = 0.5) +
            geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
            scale_color_manual(values = c("#181818", pal_biopipen()(1)), breaks = c(TRUE, FALSE)) +
            labs(x = "Sample", y = feat) +
            theme_minimal()

        vlnplot = file.path(plotsdir, paste0(slugify(feat), ".vln.png"))
        png(
            vlnplot,
            width = 800 + length(samples) * 15, height = 600, res = 100
        )
        print(vln_p)
        dev.off()

        add_report(
            list(
                src = vlnplot,
                name = feat,
                descr = paste0("Distribution of ", feat, " for each sample.")
            ),
            h1 = "Violin Plots",
            ui = "table_of_images"
        )
    }

    # Scatter plots against nCount_RNA
    log_info("- Plotting scatter plots ...")
    add_report(
        list(
            kind = "descr",
            content = paste(
                "The scatter plots for each feature against nCount_RNA. ",
                "The cells that fail the QC criteria are colored in red, and",
                "the cells that pass the QC criteria are colored in black.",
                "The cells that fail the QC criteria are filtered out in the returned Seurat object."
            )
        ),
        h1 = "Scatter Plots"
    )
    for (feat in setdiff(feats, "nCount_RNA")) {
        log_info("  For feature: {feat}, against nCount_RNA")
        scat_p <- ggplot(cell_qc_df, aes(x = nCount_RNA, y = !!sym(feat), color = .QC)) +
            geom_point() +
            scale_color_manual(values = c("#181818", pal_biopipen()(1)), breaks = c(TRUE, FALSE)) +
            labs(x = "nCount_RNA", y = feat) +
            theme_minimal()

        scatfile = file.path(plotsdir, paste0(slugify(feat), "-nCount_RNA.scatter.png"))
        png(scatfile, width = 800, height = 600, res = 100)
        print(scat_p)
        dev.off()

        add_report(
            list(
                src = scatfile,
                name = paste0(feat, " vs nCount_RNA"),
                descr = paste0("Scatter plot for ", feat, " against nCount_RNA")
            ),
            h1 = "Scatter Plots",
            ui = "table_of_images"
        )
    }

    # return the dim_df calculated from the cell_qc_df
    rbind(
        cell_qc_df %>%
            # group_by(Sample) %>%
            summarise(
                when = "Before_Cell_QC",
                nCells = dplyr::n(),
                nGenes = ngenes
            ) %>%
            ungroup(),
        cell_qc_df %>%
            filter(.QC) %>%
            # group_by(Sample) %>%
            summarise(
                when = "After_Cell_QC",
                nCells = dplyr::n(),
                nGenes = ngenes
            ) %>%
            ungroup()
    )
}

load_sample = function(sample) {
    log_info("- Loading sample: {sample} ...")
    mdata = as.data.frame(metadata)[metadata$Sample == sample, , drop=TRUE]
    path = as.character(mdata$RNAData)
    if (is.na(path) || !is.character(path) || nchar(path) == 0 || path == "NA") {
        warning(paste0("No path found for sample: ", sample))
        return (NULL)
    }

    # obj_list = list()
    if (dir.exists(path)) {
        exprs = tryCatch(
            # Read10X requires
            # - barcodes.tsv.gz
            # - genes.tsv.gz
            # - matrix.mtx.gz
            # But sometimes, they are prefixed with sample name
            # e.g.GSM4143656_SAM24345863-ln1.barcodes.tsv.gz
            { Read10X(data.dir = path) },
            error = function(e) rename_files(e, sample, path)
        )
    } else {
        exprs = Read10X_h5(path)
    }
    if ("Gene Expression" %in% names(exprs)) {
        exprs = exprs[["Gene Expression"]]
    }
    obj <- CreateSeuratObject(exprs, project=sample)
    # filter the cells that don't have any gene expressions
    # cell_exprs = colSums(obj@assays$RNA)
    # obj = subset(obj, cells = names(cell_exprs[cell_exprs > 0]))
    obj = RenameCells(obj, add.cell.id = sample)
    # Attach meta data
    for (mname in names(mdata)) {
        if (mname %in% c("RNAData", "TCRData")) { next }
        mdt = mdata[[mname]]
        if (is.factor(mdt)) { mdt = levels(mdt)[mdt] }
        obj[[mname]] = mdt
    }

    if (isTRUE(envs$cell_qc_per_sample)) {
        log_info("- Perform cell QC for sample: {sample} ...")
        obj = perform_cell_qc(obj, per_sample = TRUE)
    }

    if (isTRUE(envs$use_sct)) {
        # so that we have data and scale.data layers on RNA assay
        # useful for visualization in case some genes are not in
        # the SCT assay
        obj = NormalizeData(obj, verbose = FALSE)
        obj = FindVariableFeatures(obj, verbose = FALSE)
        obj = ScaleData(obj, verbose = FALSE)
    }
    obj
}

run_gene_qc <- function(sobj) {
    cached <- get_cached(
        list(
            cell_qc = envs$cell_qc,
            gene_qc = envs$gene_qc,
            cell_qc_per_sample = envs$cell_qc_per_sample,
            use_sct = envs$use_sct
        ),
        "GeneQC",
        cache_dir
    )
    if (!is.null(cached$data)) {
        log_info("Loading gene-QC'ed object from cache ...")
        sobj <- cached$data
    } else {
        log_info("Filtering genes ...")
        genes <- rownames(sobj)
        filtered <- FALSE
        if (!is.null(envs$gene_qc$min_cells) && envs$gene_qc$min_cells > 0) {
            genes = genes[Matrix::rowSums(sobj) >= envs$gene_qc$min_cells]
            filtered <- TRUE
        }
        excludes <- envs$gene_qc$excludes
        if (!is.null(excludes)) {
            if (length(excludes) == 1) {
                excludes <- trimws(unlist(strsplit(excludes, ",")))
            }
            for (ex in excludes) {
                genes <- genes[!grepl(ex, genes)]
            }
            filtered <- TRUE
        }
        if (filtered) {
            sobj = subset(sobj, features = genes)
        }
        cached$data <- sobj
        save_to_cache(cached, "GeneQC", cache_dir)
    }
    sobj
}

run_cell_qc <- function(sobj) {
    cached <- get_cached(
        list(cell_qc = envs$cell_qc, cell_qc_per_sample = envs$cell_qc_per_sample, use_sct = envs$use_sct),
        "CellQC",
        cache_dir
    )
    if (!is.null(cached$data)) {
        log_info("Loading cell-QC'ed object from cache ...")
        sobj <- cached$data$sobj
        cell_qc_df <<- cached$data$cell_qc_df
    } else {
        # Load data
        log_info("Reading samples individually ...")
        obj_list = lapply(samples, load_sample)

        log_info("Merging samples ...")
        sobj = Reduce(merge, obj_list)
        rm(obj_list)
        gc()

        if (!envs$cell_qc_per_sample) {
            log_info("Performing cell QC ...")
            sobj = perform_cell_qc(sobj, per_sample = FALSE)
        }

        cached$data <- list(sobj = sobj, cell_qc_df = cell_qc_df)
        save_to_cache(cached, "CellQC", cache_dir)
    }
    sobj
}

run_transformation <- function(sobj) {
    envs_cache <- envs
    envs_cache$ncores <- NULL
    envs_cache$doublet_detector <- NULL
    envs_cache$DoubletFinder <- NULL
    envs_cache$scDblFinder <- NULL
    envs_cache$IntegrateLayers <- NULL
    cached <- get_cached(envs_cache, "Transformed", cache_dir)
    if (!is.null(cached$data)) {
        log_info("Loading transformed object from cache ...")
        sobj <- cached$data
    } else {
        log_info("Performing transformation/scaling ...")
        # Not joined yet
        # sobj[["RNA"]] <- split(sobj[["RNA"]], f = sobj$Sample)
        if (envs$use_sct) {
            log_info("- Running SCTransform ...")
            SCTransformArgs <- envs$SCTransform
            # log to stdout but don't populate it to running log
            print(paste0("  SCTransform: ", format_args(SCTransformArgs)))
            log_debug("  SCTransform: {format_args(SCTransformArgs)}")
            SCTransformArgs$object <- sobj
            sobj <- do_call(SCTransform, SCTransformArgs)
            # Default is to use the SCT assay

            # Cleanup memory
            SCTransformArgs$object <- NULL
            rm(SCTransformArgs)
            gc()
        } else {
            log_info("- Running NormalizeData ...")
            NormalizeDataArgs <- envs$NormalizeData
            print(paste0("  NormalizeData: ", format_args(NormalizeDataArgs)))
            log_debug("  NormalizeData: {format_args(NormalizeDataArgs)}")
            NormalizeDataArgs$object <- sobj
            sobj <- do_call(NormalizeData, NormalizeDataArgs)

            # Cleanup memory
            NormalizeDataArgs$object <- NULL
            rm(NormalizeDataArgs)
            gc()

            log_info("- Running FindVariableFeatures ...")
            FindVariableFeaturesArgs <- envs$FindVariableFeatures
            print(paste0("  FindVariableFeatures: ", format_args(FindVariableFeaturesArgs)))
            log_debug("  FindVariableFeatures: {format_args(FindVariableFeaturesArgs)}")
            FindVariableFeaturesArgs$object <- sobj
            sobj <- do_call(FindVariableFeatures, FindVariableFeaturesArgs)

            # Cleanup memory
            FindVariableFeaturesArgs$object <- NULL
            rm(FindVariableFeaturesArgs)
            gc()

            log_info("- Running ScaleData ...")
            ScaleDataArgs <- envs$ScaleData
            print(paste0("  ScaleData: ", format_args(ScaleDataArgs)))
            log_debug("  ScaleData: {format_args(ScaleDataArgs)}")
            ScaleDataArgs$object <- sobj
            sobj <- do_call(ScaleData, ScaleDataArgs)

            # Cleanup memory
            ScaleDataArgs$object <- NULL
            rm(ScaleDataArgs)
            gc()
        }

        log_info("- Running RunPCA ...")
        RunPCAArgs <- envs$RunPCA
        RunPCAArgs$npcs <- if (is.null(RunPCAArgs$npcs)) { 50 } else { min(RunPCAArgs$npcs, ncol(sobj) - 1) }
        print(paste0("  RunPCA: ", format_args(RunPCAArgs)))
        log_debug("  RunPCA: {format_args(RunPCAArgs)}")
        RunPCAArgs$object <- sobj
        sobj <- do_call(RunPCA, RunPCAArgs)

        # Cleanup memory
        RunPCAArgs$object <- NULL
        rm(RunPCAArgs)
        gc()

        cached$data <- sobj
        save_to_cache(cached, "Transformed", cache_dir)
    }

    sobj
}

run_integration <- function(sobj) {

    envs_cache <- envs
    envs_cache$ncores <- NULL
    envs_cache$doublet_detector <- NULL
    envs_cache$DoubletFinder <- NULL
    envs_cache$scDblFinder <- NULL
    cached <- get_cached(envs_cache, "Integrated", cache_dir)

    if (!is.null(cached$data)) {
        log_info("Loading integrated/layer-joined object from cache ...")
        sobj <- cached$data
    } else {

        if (!envs$no_integration) {
            log_info("- Running IntegrateLayers (method = {envs$IntegrateLayers$method}) ...")
            IntegrateLayersArgs <- envs$IntegrateLayers
            method <- IntegrateLayersArgs$method
            if (!is.null(IntegrateLayersArgs$reference) && is.character(IntegrateLayersArgs$reference)) {
                log_info("  Using reference samples: {paste(IntegrateLayersArgs$reference, collapse = ', ')}")
                IntegrateLayersArgs$reference <- match(IntegrateLayersArgs$reference, samples)
                log_info("  Transferred to indices: {paste(IntegrateLayersArgs$reference, collapse = ', ')}")
            }
            if (method %in% c("CCA", "cca")) { method <- "CCAIntegration" } else
            if (method %in% c("RPCA", "rpca")) { method <- "RPCAIntegration" } else
            if (method %in% c("Harmony", "harmony")) { method <- "HarmonyIntegration" } else
            if (method %in% c("FastMNN", "fastmnn")) { method <- "FastMNNIntegration" } else
            if (method %in% c("scVI", "scvi")) { method <- "scVIIntegration" } else
            { stop(paste0("Unknown integration method: ", method)) }
            if (envs$use_sct && is.null(IntegrateLayersArgs$normalization.method)) {
                IntegrateLayersArgs$normalization.method <- "SCT"
            }
            IntegrateLayersArgs$method <- eval(parse(text = method))
            new_reductions <- list(
                "CCAIntegration" = "integrated.cca",
                "RPCAIntegration" = "integrated.rpca",
                "HarmonyIntegration" = "harmony",
                "FastMNNIntegration" = "integration.mnn",
                "scVIIntegration" = "integrated.scvi"
            )
            if (is.null(IntegrateLayersArgs$new.reduction)) {
                IntegrateLayersArgs$new.reduction <- new_reductions[[method]]
            }
            print(paste0("  IntegrateLayers: ", format_args(IntegrateLayersArgs)))
            log_debug("  IntegrateLayers: {format_args(IntegrateLayersArgs)}")
            IntegrateLayersArgs$object <- sobj
            sobj <- do_call(IntegrateLayers, IntegrateLayersArgs)
            # Save it for dimension reduction plots
            sobj@misc$integrated_new_reduction <- IntegrateLayersArgs$new.reduction

            # Cleanup memory
            IntegrateLayersArgs$object <- NULL
            rm(IntegrateLayersArgs)
            gc()
        }

        if (!envs$use_sct) {
            log_info("- Joining layers ...")
            sobj <- JoinLayers(sobj)
        }

        cached$data <- sobj
        save_to_cache(cached, "Integrated", cache_dir)
    }

    sobj
}