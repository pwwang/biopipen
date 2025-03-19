library(Seurat)
library(future)
library(bracer)
library(dplyr)
library(glue)
library(biopipen.utils)

metafile <- {{in.metafile | quote}}
rdsfile <- {{out.rdsfile | quote}}
joboutdir <- {{job.outdir | quote}}
envs <- {{envs | r: todot = "-", skip = 1}}

if (isTRUE(envs$cache)) { envs$cache <- joboutdir }

log <- get_logger()
reporter <- get_reporter()

set.seed(8525)
# 8TB
options(future.globals.maxSize = 8 * 1024 ^ 4)
options(future.rng.onMisuse="ignore")
options(Seurat.object.assay.version = "v5")
plan(strategy = "multicore", workers = envs$ncores)

reporter$add(
    list(
        kind = "descr",
        name = "Filters applied",
        content = paste0(
            "<p>Cell filters: ", html_escape(envs$cell_qc), "</p>",
            "<p>Gene filters: </p>",
            "<p>- Min Cells: ", envs$gene_qc$min_cells, "</p>",
            "<p>- Excludes: ", html_escape(envs$gene_qc$excludes %||% "Not set"), "</p>"
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

meta_cols = colnames(metadata)
if (!"Sample" %in% meta_cols) {
    stop("Error: Column `Sample` is not found in metafile.")
}
if (!"RNAData" %in% meta_cols) {
    stop("Error: Column `RNAData` is not found in metafile.")
}

qcdir = file.path(joboutdir, "qc")
dir.create(qcdir, showWarnings = FALSE, recursive = TRUE)

sobj <- LoadSeuratAndPerformQC(
    metadata,
    per_sample_qc = envs$cell_qc_per_sample,
    cell_qc = envs$cell_qc,
    gene_qc = envs$gene_qc,
    tmpdir = joboutdir,
    log = log,
    cache = envs$cache)

log$info("Saving dimension table ...")
dim_df <- data.frame(
    when = c("Before QC", "After QC"),
    nCells = c(nrow(sobj@misc$cell_qc_df), sum(sobj@misc$cell_qc_df$.QC)),
    nGenes = c(sobj@misc$gene_qc$before, sobj@misc$gene_qc$after)
)
write.table(dim_df, file = file.path(qcdir, "dim.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

reporter$add(
    list(
        kind = "descr",
        content = "The dimension table for the Seurat object. The table contains the number of cells and genes before and after QC. Note that the cell QC is performed before gene QC."
    ),
    list(
        kind = "table",
        data = list(path = file.path(qcdir, "dim.txt"))
    ),
    h1 = "Filters and QC",
    h2 = "Dimension table"
)

log$info("Visualizing QC metrics ...")
for (pname in names(envs$qc_plots)) {
    args <- envs$qc_plots[[pname]]
    args$kind <- args$kind %||% "cell"
    args$devpars <- args$devpars %||% list()
    args$more_formats <- args$more_formats %||% character()
    args$save_code <- args$save_code %||% FALSE
    extract_vars(args, "kind", "devpars", "more_formats", "save_code")
    if (kind == "gene") kind <- "gene_qc"
    if (kind == "cell") kind <- "cell_qc"
    args$object <- sobj
    plot_fn <- if (kind == "cell_qc") {
        gglogger::register(VizSeuratCellQC)
    } else {
        gglogger::register(VizSeuratGeneQC)
    }
    p <- do_call(plot_fn, args)
    prefix <- file.path(qcdir, paste0(slugify(pname), "_", kind))
    save_plot(p, prefix, devpars, formats = c("png", more_formats))
    if (save_code) {
        save_plotcode(p, prefix,
            setup = c("library(biopipen.utils)", "load('data.RData')", "invisible(list2env('args'))"),
            "args",
            auto_data_setup = FALSE)
    }
    reporter$add(
        reporter$image(prefix, more_formats, save_code, kind = "image"),
        h1 = "Filters and QC",
        h2 = html_escape(pname)
    )
}

sobj <- RunSeuratTransformation(
    sobj,
    use_sct = envs$use_sct,
    SCTransformArgs = envs$SCTransform,
    NormalizeDataArgs = envs$NormalizeData,
    FindVariableFeaturesArgs = envs$FindVariableFeatures,
    ScaleDataArgs = envs$ScaleData,
    RunPCAArgs = envs$RunPCA,
    log = log,
    cache = envs$cache
)
sobj <- RunSeuratIntegration(
    sobj,
    no_integration = envs$no_integration,
    IntegrateLayersArgs = envs$IntegrateLayers,
    log = log,
    cache = envs$cache
)

# This is the last step, doesn't need to be cached
if (!identical(envs$doublet_detector, "none")) {
    dbldir <- file.path(joboutdir, "doublets")
    dir.create(dbldir, showWarnings = FALSE, recursive = TRUE)

    sobj <- RunSeuratDoubletDetection(
        sobj,
        tool = envs$doublet_detector,
        DoubletFinderArgs = envs$DoubletFinder,
        scDblFinderArgs = envs$scDblFinder,
        filter = FALSE,
        log = log,
        cache = envs$cache
    )

    log$info("Visualizing doublet detection results ...")
    if (identical(tolower(envs$doublet_detector), "doubletfinder")) {
        p <- VizSeuratDoublets(sobj, plot_type = "pK", x_text_angle = 90)
        save_plot(
            p, file.path(dbldir, "doubletfinder_pk"),
            devpars = list(res = 100, width = 800, height = 600),
            formats = "png")
        reporter$add(
            list(
                kind = "descr",
                content = paste(
                    "The pK plot from DoubletFinder to select the optimal pK value.",
                    "See more at https://github.com/chris-mcginnis-ucsf/DoubletFinder"
                )
            ),
            list(
                kind = "image",
                src = file.path(dbldir, "doubletfinder_pk.png")
            ),
            h1 = glue("Doublet detection using {envs$doublet_detector}"),
            h2 = "BC metric vs pK"
        )
    }

    for (pt in c("dim", "pie")) {
        p <- VizSeuratDoublets(sobj, plot_type = pt)
        save_plot(p, file.path(dbldir, paste0("doublets_", pt)), formats = "png")

        reporter$add(
            list(
                src = file.path(dbldir, paste0("doublets_", pt, ".png")),
                descr = ifelse(pt == "dim", "Dimention Reduction Plot", "Pie Chart")
            ),
            h1 = glue("Doublet detection using {envs$doublet_detector}"),
            h2 = "Doublets distribution",
            ui = "table_of_images"
        )
    }

    sobj <- subset(sobj, subset = !!sym(paste0(sobj@misc$doublets$tool, "_DropletType")) != "doublet")
}

log$info("Saving QC'ed seurat object ...")
reporter$save(joboutdir)
saveRDS(sobj, rdsfile)
