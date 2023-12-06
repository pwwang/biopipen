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
envs = {{envs | r}}

set.seed(8525)
options(future.globals.maxSize = 80000 * 1024^2)
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
    if (is.na(path) || !is.character(path) || nchar(path) == 0) {
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
    obj = CreateSeuratObject(counts=exprs, project=sample)
    # filter the cells that don't have any gene expressions
    cell_exprs = colSums(obj@assays$RNA)
    obj = subset(obj, cells = names(cell_exprs[cell_exprs > 0]))
    # obj = SCTransform(object=obj, return.only.var.genes=FALSE, verbose=FALSE)
    obj = RenameCells(obj, add.cell.id = sample)
    # Attach meta data
    for (mname in names(mdata)) {
        if (mname %in% c("RNAData", "TCRData")) { next }
        mdt = mdata[[mname]]
        if (is.factor(mdt)) { mdt = levels(mdt)[mdt] }
        obj[[mname]] = mdt
    }
    # obj_list[[sample]] = obj

    # obj_list
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
    if ("min_cells" %in% names(envs$gene_qc)) {
        genes = rownames(sobj)[Matrix::rowSums(sobj) >= envs$gene_qc$min_cells]
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


log_info("Saving filtered seurat object ...")
saveRDS(sobj, rdsfile)

save_report(joboutdir)
