library(rlang)
library(dplyr)
library(plotthis)
library(biopipen.utils)

indirs <- {{in.indirs | r}}
outdir <- {{out.outdir | r}}
joboutdir <- {{job.outdir | r}}
metrics_glob <- {{envs.metrics_glob | r}}

logger <- get_logger()
reporter <- get_reporter()

if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

logger$info("Looking for metrics files matching '{metrics_glob}' ...")
metrics_files <- data.frame(
    Sample = character(0),
    File = character(0),
    Rel = character(0)
)
for (indir in indirs) {
    sample <- basename(indir)
    files <- list.files(
        indir,
        pattern = glob2rx(metrics_glob),
        recursive = TRUE,
        full.names = TRUE
    )
    if (length(files) == 0) {
        logger$warn("- No metrics file matching '{metrics_glob}' in: {indir}")
        next
    }
    files <- sort(files)
    logger$info("- Found {length(files)} metrics file(s) in: {indir}")
    metrics_files <- rbind(
        metrics_files,
        data.frame(
            Sample = sample,
            File = files,
            Rel = sub(paste0("^", indir, "/*"), "", files)
        )
    )
}

if (nrow(metrics_files) == 0) {
    stop(paste0(
        "No metrics file matching '", metrics_glob, "' was found in any of ",
        "the input directories. Make sure `envs.metrics_glob` is correct ",
        "for your version of split-pipe. Input directories: ",
        paste0(indirs, collapse = ", ")
    ))
}

# Keep a copy of the raw metrics files
rawdir <- file.path(outdir, "raw")
logger$info("Copying raw metrics files to: {rawdir}")
for (i in seq_len(nrow(metrics_files))) {
    dest <- file.path(rawdir, metrics_files$Sample[i], metrics_files$Rel[i])
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(metrics_files$File[i], dest, overwrite = TRUE)
}

logger$info("Reading and merging metrics ...")
metrics <- NULL
for (i in seq_len(nrow(metrics_files))) {
    sample <- metrics_files$Sample[i]
    file <- metrics_files$File[i]
    logger$debug("- Reading metrics: {file}")
    metric <- read.csv(file, header = TRUE, row.names = NULL, check.names = FALSE)
    metric <- cbind(
        data.frame(Sample = sample, File = metrics_files$Rel[i]),
        metric
    )
    if (!is.null(metrics)) {
        missing_cols <- setdiff(colnames(metrics), colnames(metric))
        if (length(missing_cols) > 0) {
            logger$warn('- Missing columns: {paste0(missing_cols, collapse = ", ")} in: {file}')
            metric[missing_cols] <- NA
        }
        missing_cols <- setdiff(colnames(metric), colnames(metrics))
        if (length(missing_cols) > 0) {
            logger$warn('- Missing columns: {paste0(missing_cols, collapse = ", ")} in metrics before: {file}')
            metrics[missing_cols] <- NA
        }
    }
    metrics <- rbind(metrics, metric)
}

# Save metrics
write.table(
    metrics,
    file.path(outdir, "metrics.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

reporter$add(
    list(
        kind = "descr",
        content = paste0(
            "Metrics of all sublibraries merged into one table ",
            "(`raw/` keeps a copy of each metrics file). ",
            "One section per metric with its description and a bar plot ",
            "follows the table."
        )
    ),
    list(kind = "table", src = file.path(outdir, "metrics.txt")),
    h1 = "Metrics of all sublibraries"
)

# Column descriptions for the merged metrics. split-pipe's own documentation is
# behind a login, so these are written from its outputs: the *_Q30 and *_fraction
# columns are fractions in [0, 1], TSCP is split-pipe's molecule count, and the
# test_* columns are its run-level QC statistics.
METRIC_DESCR <- list(
    `Sample` = "The input directory this row came from (the sublibrary or the combined run).",
    `File` = "The metrics file the row was read from, relative to that directory.",
    `sample` = "split-pipe's sample name (usually `all-sample` for an all-well run).",
    `sample_well_count` = "Number of wells defined for the sample.",
    `number_of_cells` = "Cells called for the sample.",
    `test_median_tscp_per_cell` = "Median TSCPs per cell.",
    `test_median_tscp_at50` = "Median TSCPs per cell at half depth (NA when split-pipe does not compute it).",
    `test_median_genes_per_cell` = "Median genes per cell.",
    `mean_reads_per_cell` = "Number of input reads divided by the number of cells.",
    `number_of_reads` = "Total number of input reads.",
    `number_of_tscp` = "Total number of TSCPs (molecules) counted.",
    `sequencing_saturation` = "Fraction of reads originating from an already-observed molecule (UMI): 0 means every read was a new molecule.",
    `valid_barcode_fraction` = "Fraction of reads whose barcodes match the whitelist after correction.",
    `bc1_Q30` = "Fraction of round-1 barcode bases with Q-score >= 30.",
    `bc2_Q30` = "Fraction of round-2 barcode bases with Q-score >= 30.",
    `bc3_Q30` = "Fraction of round-3 barcode bases with Q-score >= 30.",
    `cDNA_Q30` = "Fraction of cDNA (read 1) bases with Q-score >= 30.",
    `polyN_Q30` = "Fraction of the R2 poly-N (UMI block) bases with Q-score >= 30.",
    `sample_tmap_fraction` = "Fraction of the sample's reads that map to the transcriptome.",
    `transcriptome_map_fraction` = "Fraction of reads that map to the transcriptome.",
    `tso_fraction_in_read1` = "Fraction of read-1 reads containing the TSO sequence (library-quality indicator; 0 for a clean library).",
    `test_fraction_reads_in_cells` = "Fraction of reads that fall in called cells.",
    `test_fraction_tscp_in_cells` = "Fraction of TSCPs that fall in called cells.",
    `test_fraction_exonic` = "Fraction of aligned reads that are exonic.",
    `cell_tscp_cutoff` = "The TSCP cutoff used to call cells (from `dge_tscp_cut_given`, or estimated)."
)

logger$info("Plotting metrics ...")
for (metric in colnames(metrics)) {
    if (metric %in% c("Sample", "File")) next
    if (!is.numeric(metrics[[metric]])) {
        logger$info("- {metric}: not numeric, skipping plot")
        next
    }
    logger$info("- {metric}")

    reporter$add(
        list(
            kind = "descr",
            content = METRIC_DESCR[[metric]] %||% paste0("Metric: ", metric)
        ),
        h1 = metric
    )

    # barplot over the input directories
    p <- BarPlot(metrics, x = "Sample", y = metric, x_text_angle = 90)
    figfile <- file.path(outdir, paste0(slugify(metric), ".barplot.png"))
    png(figfile, height = 600, res = 100, width = max(nrow(metrics) * 30 + 200, 400))
    print(p)
    dev.off()

    reporter$add(
        list(src = figfile, name = "By sample"),
        ui = "table_of_images",
        h1 = metric
    )
}

reporter$save(joboutdir)
