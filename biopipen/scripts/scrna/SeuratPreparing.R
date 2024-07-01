source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/caching.R")

library(Seurat)
library(future)
library(bracer)
library(ggplot2)
library(dplyr)
# library(tidyseurat)

metafile <- {{in.metafile | quote}}
rdsfile <- {{out.rdsfile | quote}}
joboutdir <- {{job.outdir | quote}}
envs <- {{envs | r: todot = "-", skip = 1}}

if (isTRUE(envs$cache)) { envs$cache <- joboutdir }
if (length(envs$cache) > 1) {
    log_warn("Multiple cache directories (envs.cache) detected, using the first one.")
    envs$cache <- envs$cache[1]
}

set.seed(8525)
# 8TB
options(future.globals.maxSize = 8 * 1024 ^ 4)
options(future.rng.onMisuse="ignore")
options(Seurat.object.assay.version = "v5")
plan(strategy = "multicore", workers = envs$ncores)

.stringify_list <- function(x) {
    paste(sapply(names(x), function(n) paste(n, x[[n]], sep = " = ") ), collapse = "; ")
}

add_report(
    list(
        kind = "descr",
        name = "Filters applied",
        content = paste0(
            "<p>Cell filters: ", html_escape(envs$cell_qc), "</p>",
            "<p>Gene filters: ", html_escape(.stringify_list(envs$gene_qc)), "</p>"
        )
    ),
    h1 = "Filters and QC"
)

metadata <- read.table(
    metafile,
    header = TRUE,
    row.names = NULL,
    sep = "\t",
    check.names = FALSE
)

cache_sig <- capture.output(str(metadata))
dig_sig <- digest::digest(cache_sig, algo = "md5")
dig_sig <- substr(dig_sig, 1, 8)
cache_dir <- NULL
if (is.character(envs$cache)) {
    cache_dir <- file.path(envs$cache, paste0(dig_sig, ".seuratpreparing_cache"))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(cache_sig, file.path(cache_dir, "signature.txt"))
}

meta_cols = colnames(metadata)
if (!"Sample" %in% meta_cols) {
    stop("Error: Column `Sample` is not found in metafile.")
}
if (!"RNAData" %in% meta_cols) {
    stop("Error: Column `RNAData` is not found in metafile.")
}

samples = as.character(metadata$Sample)

# used for plotting
cell_qc_df = NULL

plotsdir = file.path(joboutdir, "plots")
dir.create(plotsdir, showWarnings = FALSE, recursive = TRUE)

# features for cell QC
feats = c(
    "nFeature_RNA", "nCount_RNA",
    "percent.mt", "percent.ribo", "percent.hb", "percent.plat"
)

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
    sobj$percent.mt <- PercentageFeatureSet(sobj, pattern = "^MT-")
    sobj$percent.ribo <- PercentageFeatureSet(sobj, pattern = "^RP[SL]")
    sobj$percent.hb <- PercentageFeatureSet(sobj, pattern = "^HB[^(P)]")
    sobj$percent.plat <- PercentageFeatureSet(sobj, pattern = "PECAM1|PF4")

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
        obj = perform_cell_qc(obj, TRUE)
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

cached <- get_cached(
    list(cell_qc = envs$cell_qc, cell_qc_per_sample = envs$cell_qc_per_sample, use_sct = envs$use_sct),
    "CellQC",
    cache_dir
)
if (!is.null(cached$data)) {
    log_info("Loading cell-QC'ed object from cache ...")
    sobj <- cached$data$sobj
    cell_qc_df <- cached$data$cell_qc_df
    cached$data$sobj <- NULL
    cached$data$cell_qc_df <- NULL
    cached$data <- NULL
    rm(cached)
    gc()
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
        sobj = perform_cell_qc(sobj)
    }

    cached$data = list(sobj = sobj, cell_qc_df = cell_qc_df)
    save_to_cache(cached, "CellQC", cache_dir)
}

# plot and report the QC
log_info("Plotting and reporting QC ...")
dim_df = report_cell_qc(nrow(sobj))

if (is.list(envs$gene_qc)) {
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
        cached$data <- NULL
        rm(cached)
        gc()
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
}
dim_df = rbind(
    dim_df,
    data.frame(
        when = "After_Gene_QC",
        nCells = ncol(sobj),
        nGenes = nrow(sobj)
    )
)

log_info("Saving dimension table ...")
write.table(dim_df, file = file.path(plotsdir, "dim.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

add_report(
    list(
        kind = "descr",
        content = paste(
            "The dimension table for the Seurat object. The table contains the number of cells and genes before and after QC."
        )
    ),
    list(
        kind = "table",
        data = list(path = file.path(plotsdir, "dim.txt"))
    ),
    h1 = "Filters and QC"
)

.formatArgs <- function(args) {
    paste(capture.output(str(args)), collapse = ", ")
}

envs_cache <- envs
envs_cache$ncores <- NULL
envs_cache$DoubletFinder <- NULL
envs_cache$IntegrateLayers <- NULL
cached <- get_cached(envs_cache, "Transformed", cache_dir)
if (!is.null(cached$data)) {
    log_info("Loading transformed object from cache ...")
    sobj <- cached$data
    cached$data <- NULL
    rm(cached)
    gc()
} else {
    log_info("Performing transformation/scaling ...")
    # Not joined yet
    # sobj[["RNA"]] <- split(sobj[["RNA"]], f = sobj$Sample)
    if (envs$use_sct) {
        log_info("- Running SCTransform ...")
        SCTransformArgs <- envs$SCTransform
        # log to stdout but don't populate it to running log
        print(paste0("  SCTransform: ", .formatArgs(SCTransformArgs)))
        log_debug("  SCTransform: {.formatArgs(SCTransformArgs)}")
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
        print(paste0("  NormalizeData: ", .formatArgs(NormalizeDataArgs)))
        log_debug("  NormalizeData: {.formatArgs(NormalizeDataArgs)}")
        NormalizeDataArgs$object <- sobj
        sobj <- do_call(NormalizeData, NormalizeDataArgs)

        # Cleanup memory
        NormalizeDataArgs$object <- NULL
        rm(NormalizeDataArgs)
        gc()

        log_info("- Running FindVariableFeatures ...")
        FindVariableFeaturesArgs <- envs$FindVariableFeatures
        print(paste0("  FindVariableFeatures: ", .formatArgs(FindVariableFeaturesArgs)))
        log_debug("  FindVariableFeatures: {.formatArgs(FindVariableFeaturesArgs)}")
        FindVariableFeaturesArgs$object <- sobj
        sobj <- do_call(FindVariableFeatures, FindVariableFeaturesArgs)

        # Cleanup memory
        FindVariableFeaturesArgs$object <- NULL
        rm(FindVariableFeaturesArgs)
        gc()

        log_info("- Running ScaleData ...")
        ScaleDataArgs <- envs$ScaleData
        print(paste0("  ScaleData: ", .formatArgs(ScaleDataArgs)))
        log_debug("  ScaleData: {.formatArgs(ScaleDataArgs)}")
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
    print(paste0("  RunPCA: ", .formatArgs(RunPCAArgs)))
    log_debug("  RunPCA: {.formatArgs(RunPCAArgs)}")
    RunPCAArgs$object <- sobj
    sobj <- do_call(RunPCA, RunPCAArgs)

    # Cleanup memory
    RunPCAArgs$object <- NULL
    rm(RunPCAArgs)
    gc()

    cached$data <- sobj
    save_to_cache(cached, "Transformed", cache_dir)
}

envs_cache <- envs
envs_cache$ncores <- NULL
envs_cache$DoubletFinder <- NULL
cached <- get_cached(envs_cache, "Integrated", cache_dir)

if (!is.null(cached$data)) {
    log_info("Loading integrated/layer-joined object from cache ...")
    sobj <- cached$data
    cached$data <- NULL
    rm(cached)
    gc()

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
        print(paste0("  IntegrateLayers: ", .formatArgs(IntegrateLayersArgs)))
        log_debug("  IntegrateLayers: {.formatArgs(IntegrateLayersArgs)}")
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


# This is the last step, doesn't need to be cached
if (!is.null(envs$DoubletFinder) && is.list(envs$DoubletFinder) && envs$DoubletFinder$PCs > 0) {
    library(DoubletFinder)

    log_info("Running DoubletFinder ...")
    log_info("- Preparing Seurat object ...")

    if (is.null(envs$DoubletFinder$ncores)) {
        envs$DoubletFinder$ncores <- envs$ncores
    }

    # More controls from envs?
    sobj <- FindNeighbors(sobj, dims = 1:envs$DoubletFinder$PCs)
    sobj <- FindClusters(sobj)

    log_info("- pK Indentification ...")
    sweep.res.list <- paramSweep(
        sobj,
        PCs = 1:envs$DoubletFinder$PCs,
        sct = envs$use_sct,
        num.cores = envs$DoubletFinder$ncores
    )
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)

    bcmvn$Selected <- bcmvn$pK == bcmvn$pK[which.max(bcmvn$BCmetric)[1]]
    plot <- ggplot(bcmvn, aes(x = pK, y = BCmetric, color = Selected)) +
        geom_point() +
        # rotate x axis labels
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(plot, filename = file.path(plotsdir, "pK_BCmetric.png"))

    pK <- bcmvn$pK[which.max(bcmvn$BCmetric)[1]]
    pK <- as.numeric(as.character(pK))
    pN <- envs$DoubletFinder$pN
    log_info("- Homotypic Doublet Proportion Estimate ...")
    homotypic.prop <- modelHomotypic(Idents(sobj))
    nExp_poi <- round(nrow(sobj@meta.data) * envs$DoubletFinder$doublets)
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

    log_info("- Running DoubletFinder ...")
    sobj <- doubletFinder(
        sobj,
        PCs = 1:envs$DoubletFinder$PCs,
        pN = pN,
        pK = pK,
        nExp = nExp_poi.adj,
        reuse.pANN = FALSE,
        sct = envs$use_sct
    )
    pANN_col <- paste0("pANN_", pN, "_", pK)
    pANN_col <- colnames(sobj@meta.data)[grepl(pANN_col, colnames(sobj@meta.data))]
    DF_col <- paste0("DF.classifications_", pN, "_", pK)
    DF_col <- colnames(sobj@meta.data)[grepl(DF_col, colnames(sobj@meta.data))]
    doublets <- as.data.frame(
        cbind(
            colnames(sobj),
            sobj@meta.data[, pANN_col],
            sobj@meta.data[, DF_col]
        )
    )
    colnames(doublets) <-  c("Barcode","DoubletFinder_score","DoubletFinder_DropletType")
    write.table(
        doublets,
        file.path(joboutdir, "DoubletFinder_doublets_singlets.txt"),
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )

    summary <- as.data.frame(table(doublets$DoubletFinder_DropletType))
    colnames(summary) <- c("Classification", "Droplet_N")
    write.table(
        summary,
        file.path(joboutdir, "DoubletFinder_summary.txt"),
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )

    # Do a dimplot
    log_info("- Plotting dimension reduction ...")
    dimp <- DimPlot(
        sobj, group.by = DF_col, order = "Doublet",
        cols = c("#333333", "#FF3333"), pt.size = 0.8, alpha = 0.5)
    ggsave(dimp, filename = file.path(plotsdir, "DoubletFinder_dimplot.png"))

    log_info("- Filtering doublets ...")
    sobj <- subset(sobj, cells = doublets$Barcode[doublets$DoubletFinder_DropletType == "Singlet"])

    add_report(
        list(
            kind = "descr",
            content = "The table contains the number of cells classified as singlets and doublets."
        ),
        list(
            kind = "table",
            data = list(path = file.path(joboutdir, "DoubletFinder_summary.txt"))
        ),
        h1 = "DoubletFinder Results",
        h2 = "The DoubletFinder Summary"
    )
    add_report(
        list(
            name = "pK vs BCmetric",
            src = file.path(plotsdir, "pK_BCmetric.png")
        ),
        list(
            name = "Dimension Reduction Plot",
            src = file.path(plotsdir, "DoubletFinder_dimplot.png")
        ),
        ui = "table_of_images",
        h1 = "DoubletFinder Results",
        h2 = "Plots"
    )
}

log_info("Saving QC'ed seurat object ...")
saveRDS(sobj, rdsfile)

save_report(joboutdir)
