source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(future)
library(bracer)
library(ggplot2)
library(tidyseurat)
library(slugify)

metafile = {{in.metafile | quote}}
rdsfile = {{out.rdsfile | quote}}
joboutdir = {{job.outdir | quote}}
envs = {{envs | r: todot = "-", skip = 1}}

set.seed(8525)
options(future.globals.maxSize = 80000 * 1024^2)
options(Seurat.object.assay.version = "v5")
plan(strategy = "multicore", workers = envs$ncores)

add_report(
    list(
        kind = "descr",
        name = "Filters applied",
        content = paste0(
            "<p>Cell filters: ", html_escape(envs$cell_qc), "</p>",
            "<p>Gene filters: ", html_escape(envs$gene_qc), "</p>"
        )
    ),
    h1 = "Filters and QC"
)

metadata = read.table(
    metafile,
    header = TRUE,
    row.names = NULL,
    sep = "\t",
    check.names = FALSE
)

meta_cols = colnames(metadata)
if (!"Sample" %in% meta_cols) {
    stop("Error: Column `Sample` is not found in metafile.")
}
if (!"RNAData" %in% meta_cols) {
    stop("Error: Column `RNAData` is not found in metafile.")
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
    cell_exprs = colSums(obj@assays$RNA)
    obj = subset(obj, cells = names(cell_exprs[cell_exprs > 0]))
    obj = RenameCells(obj, add.cell.id = sample)
    # Attach meta data
    for (mname in names(mdata)) {
        if (mname %in% c("RNAData", "TCRData")) { next }
        mdt = mdata[[mname]]
        if (is.factor(mdt)) { mdt = levels(mdt)[mdt] }
        obj[[mname]] = mdt
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

# Load data
samples = as.character(metadata$Sample)

log_info("Reading samples individually ...")
obj_list = lapply(samples, load_sample)

log_info("Merging samples ...")
if (length(obj_list) >= 2) {
    y = c()
    for (i in 2:length(obj_list)) y = c(y, obj_list[[i]])
    sobj = merge(obj_list[[1]], y)
} else {
    sobj = obj_list[[1]]
}

log_info("Adding metadata for QC ...")
sobj$percent.mt = PercentageFeatureSet(sobj, pattern = "^MT-")
sobj$percent.ribo = PercentageFeatureSet(sobj, pattern = "^RP[SL]")
sobj$percent.hb = PercentageFeatureSet(sobj, pattern = "^HB[^(P)]")
sobj$percent.plat = PercentageFeatureSet(sobj, pattern = "PECAM1|PF4")

dim_df = data.frame(When = "Before_QC", nCells = ncol(sobj), nGenes = nrow(sobj))

if (is.null(envs$cell_qc) || length(envs$cell_qc) == 0) {
    log_warn("No cell QC criteria is provided. All cells will be kept.")
    envs$cell_qc = "TRUE"
}

sobj = sobj %>% mutate(.QC = !!rlang::parse_expr(envs$cell_qc))
feats = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat")
plotsdir = file.path(joboutdir, "plots")
dir.create(plotsdir, showWarnings = FALSE)

# Violin plots
log_info("Plotting violin plots ...")
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
    log_info("- For feature: {feat}")
    vln_p = VlnPlot(
        sobj,
        cols = rep("white", length(samples)),
        group.by = "Sample",
        features = feat,
        pt.size = 0) + NoLegend()
    vln_p$data$.QC = sobj@meta.data$.QC
    vln_p = vln_p + geom_jitter(
            aes(color = .QC),
            data = vln_p$data,
            position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9)
        ) + scale_color_manual(values = c("#181818", pal_biopipen()(1)), breaks = c(TRUE, FALSE))

    vlnplot = file.path(plotsdir, paste0(slugify(feat, tolower = FALSE), ".vln.png"))
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
log_info("Plotting scatter plots ...")
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
    log_info("- For feature: {feat}, against nCount_RNA")
    scat_p = FeatureScatter(
        sobj,
        feature1 = "nCount_RNA",
        feature2 = feat,
        group.by = ".QC"
    ) +
    NoLegend() +
    scale_color_manual(values = c("#181818", pal_biopipen()(1)), breaks = c(TRUE, FALSE))

    scatfile = file.path(plotsdir, paste0(slugify(feat, tolower = FALSE), "-nCount_RNA.scatter.png"))
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

# Do the filtering
log_info("Filtering cells using QC criteria ...")
sobj = sobj %>% filter(.QC)
sobj$.QC = NULL

log_info("Filtering genes ...")
if (is.list(envs$gene_qc)) {
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
}
dim_df = rbind(
    dim_df,
    data.frame(
        When = "After_Gene_QC",
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

log_info("Perform integration ...")
# Not joined yet
# sobj[["RNA"]] <- split(sobj[["RNA"]], f = sobj$Sample)
if (envs$use_sct) {
    log_info("- Running SCTransform ...")
    SCTransformArgs <- envs$SCTransform
    SCTransformArgs$object <- sobj
    sobj <- do_call(SCTransform, SCTransformArgs)
    # Default is to use the SCT assay
} else {
    log_info("- Running NormalizeData ...")
    NormalizeDataArgs <- envs$NormalizeData
    NormalizeDataArgs$object <- sobj
    sobj <- do_call(NormalizeData, NormalizeDataArgs)

    log_info("- Running FindVariableFeatures ...")
    FindVariableFeaturesArgs <- envs$FindVariableFeatures
    FindVariableFeaturesArgs$object <- sobj
    sobj <- do_call(FindVariableFeatures, FindVariableFeaturesArgs)

    log_info("- Running ScaleData ...")
    ScaleDataArgs <- envs$ScaleData
    ScaleDataArgs$object <- sobj
    sobj <- do_call(ScaleData, ScaleDataArgs)
}

log_info("- Running RunPCA ...")
RunPCAArgs <- envs$RunPCA
RunPCAArgs$npcs <- if (is.null(RunPCAArgs$npcs)) { 50 } else { min(RunPCAArgs$npcs, ncol(sobj) - 1) }
RunPCAArgs$object <- sobj
sobj <- do_call(RunPCA, RunPCAArgs)

if (!envs$no_integration) {
    log_info("- Running IntegrateLayers ...")
    IntegrateLayersArgs <- envs$IntegrateLayers
    IntegrateLayersArgs$object <- sobj
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
    sobj <- do_call(IntegrateLayers, IntegrateLayersArgs)
    # Save it for dimension reduction plots
    sobj@misc$integrated_new_reduction <- IntegrateLayersArgs$new.reduction
}

if (!envs$use_sct) {
    log_info("- Joining layers ...")
    sobj <- JoinLayers(sobj)
}

log_info("Saving filtered seurat object ...")
saveRDS(sobj, rdsfile)

save_report(joboutdir)
