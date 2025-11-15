library(Seurat)
library(future)
library(bracer)
library(dplyr)
library(glue)
library(biopipen.utils)

metafile <- {{in.metafile | r}}
outfile <- {{out.outfile | r}}
joboutdir <- {{job.outdir | r}}
envs <- {{envs | r: todot = "-", skip = 1}}

if (isTRUE(envs$cache)) { envs$cache <- joboutdir }

log <- get_logger()
reporter <- get_reporter()

set.seed(8525)
# 8TB
options(future.globals.maxSize = Inf)
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
            "<p>- Excludes: ",
            ifelse(is.null(envs$gene_qc$excludes), "Not set", paste(envs$gene_qc$excludes, collapse = ", ")),
            "</p>"
        )
    ),
    h1 = "Filters and QC"
)

metadata <- tryCatch({
    log$debug("Trying to read Seurat object from metafile ...")
    read_obj(metafile)
}, error = function(e) {
    log$debug("Failed to read Seurat object from metafile: {e$message}")
    log$debug("Reading metafile as a table (sample info) ...")
    read.table(
        metafile,
        header = TRUE,
        row.names = NULL,
        sep = "\t",
        check.names = FALSE
    )
})
is_seurat <- inherits(metadata, "Seurat")

meta_cols <- if (is_seurat) colnames(metadata@meta.data) else colnames(metadata)
if (!"Sample" %in% meta_cols) {
    stop("Error: Column `Sample` is not found in ", ifelse(is_seurat, "Seurat object's meta.data.", "metafile."))
}
if (!"RNAData" %in% meta_cols && !is_seurat) {
    stop("Error: Column `RNAData` is not found in metafile.")
}

qcdir = file.path(joboutdir, "qc")
dir.create(qcdir, showWarnings = FALSE, recursive = TRUE)

sobj <- LoadSeuratAndPerformQC(
    metadata,
    min_cells = envs$min_cells,
    min_features = envs$min_features,
    cell_qc = envs$cell_qc,
    gene_qc = envs$gene_qc,
    tmpdir = joboutdir,
    log = log,
    cache = envs$cache)

log$info("Saving and visualizing QC results ...")
cell_qc_df <- VizSeuratCellQC(sobj, plot_type = "table")
write.table(cell_qc_df, file = file.path(qcdir, "cell_qc.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

reporter$add(
    list(
        name = "Cell QC metrics",
        contents = list(
            list(
                kind = "descr",
                content = paste0(
                    "The table below show the number of cells in each sample that failed and passed the QC filters. ",
                    "The last row shows the total number of cells that failed and passed the QC filters across all samples. "
                )
            ),
            list(kind = "table", src = file.path(qcdir, "cell_qc.txt"))
        )
    ),
    h1 = "Filters and QC",
    h2 = "Cell-level Quality Control",
    ui = "tabs"
)

gene_qc_df <- VizSeuratGeneQC(sobj, plot_type = "table")
write.table(gene_qc_df, file = file.path(qcdir, "gene_qc.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

reporter$add(
    list(
        name = "Gene QC metrics",
        contents = list(
            list(
                kind = "descr",
                content = paste0(
                    "The table below show the number of genes in each sample that failed and passed the QC filters. ",
                    "The last row shows the final number of genes that failed and passed the QC filters across all samples. ",
                    "Any gene that failed the QC filters will be excluded in the merged Seurat object."
                )
            ),
            list(kind = "table", src = file.path(qcdir, "gene_qc.txt")),
            list(kind = "list", items = list(paste0(
                "We may still end up with features slightly less than the final passed ones. ",
                "For example, when SCTransform is used, the number of features may be less than the number of genes that passed the QC filters. ",
                "This is because SCTransform selects the top N features based on variance. "
            )))
        )
    ),
    h1 = "Filters and QC",
    h2 = "Gene-level Quality Control",
    ui = "tabs"
)

for (pname in names(envs$qc_plots)) {
    if (is.null(envs$qc_plots[[pname]])) next
    log$info("- {pname} ...")
    args <- envs$qc_plots[[pname]]
    args$kind <- args$kind %||% "cell"
    args$devpars <- args$devpars %||% list()
    args$more_formats <- args$more_formats %||% character()
    args$save_code <- args$save_code %||% FALSE
    args$descr <- args$descr %||% pname
    extract_vars(args, "kind", "devpars", "more_formats", "save_code", "descr")
    if (kind == "gene") kind <- "gene_qc"
    if (kind == "cell") kind <- "cell_qc"
    args$object <- sobj
    plot_fn <- if (kind == "cell_qc") {
        gglogger::register(VizSeuratCellQC)
    } else {
        gglogger::register(VizSeuratGeneQC)
    }
    p <- do_call(plot_fn, args)
    prefix <- file.path(qcdir, paste0(slugify(pname), ".", kind))
    save_plot(p, prefix, devpars, formats = c("png", more_formats))
    if (save_code) {
        save_plotcode(p, prefix,
            setup = c("library(biopipen.utils)", "load('data.RData')", "invisible(list2env(args, envir = .GlobalEnv))"),
            "args",
            auto_data_setup = FALSE)
    }
    reporter$add(
        list(
            name = pname,
            contents = list(
                list(kind = "descr", content = descr),
                reporter$image(prefix, more_formats, save_code, kind = "image")
            )
        ),
        h1 = "Filters and QC",
        h2 = ifelse(kind == "cell_qc", "Cell-level Quality Control", "Gene-level Quality Control"),
        ui = "tabs"
    )
}

log$info("Filtering with QC criteria ...")
sobj <- FinishSeuratQC(sobj)

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

if (!is.null(envs$mutaters) && length(envs$mutaters) > 0) {
    log$info("Mutating metadata ...")
    sobj@meta.data <- sobj@meta.data %>%
        mutate(!!!lapply(envs$mutaters, rlang::parse_expr))
}

log$info("Saving QC'ed seurat object ...")
reporter$save(joboutdir)
save_obj(sobj, outfile)
