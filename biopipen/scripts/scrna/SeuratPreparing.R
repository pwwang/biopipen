{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "caching.R" | source_r }}

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

{{ biopipen_dir | joinpaths: "scripts", "scrna", "SeuratPreparing-common.R" | source_r }}

add_report(
    list(
        kind = "descr",
        name = "Filters applied",
        content = paste0(
            "<p>Cell filters: ", html_escape(envs$cell_qc), "</p>",
            "<p>Gene filters: ", html_escape(stringify_list(envs$gene_qc)), "</p>"
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

sobj <- run_cell_qc(sobj)

# plot and report the QC
log_info("Plotting and reporting QC ...")
dim_df = report_cell_qc(nrow(sobj))

if (is.list(envs$gene_qc)) {
    sobj <- run_gene_qc(sobj)
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

sobj <- run_transformation(sobj)
sobj <- run_integration(sobj)

# This is the last step, doesn't need to be cached
if (!is.null(envs$doublet_detector) && envs$doublet_detector != "none") {
    {{* biopipen_dir | joinpaths: "scripts", "scrna", "SeuratPreparing-doublet_detection.R" | source_r }}

    detector <- tolower(envs$doublet_detector)
    if (detector == "doubletfinder") detector <- "DoubletFinder"
    if (detector == "scdblfinder") detector <- "scDblFinder"
    dd <- run_dd(detector)
    save_dd(dd, detector)
    sobj <- add_dd_to_seurat(sobj, dd)
    plot_dd(sobj, dd, detector)
    sobj <- filter_dd(sobj, dd, detector)
    report_dd(detector)
}


log_info("Saving QC'ed seurat object ...")
saveRDS(sobj, rdsfile)

save_report(joboutdir)
