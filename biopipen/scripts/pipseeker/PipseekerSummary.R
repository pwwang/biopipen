library(rlang)
library(dplyr)
library(plotthis)
library(biopipen.utils)

indirs <- {{in.indirs | r}}
outdir <- {{out.outdir | r}}
joboutdir <- {{job.outdir | r}}
group <- {{envs.group | r}}
sensitivity <- {{envs.sensitivity | r}}

logger <- get_logger()
reporter <- get_reporter()

if (is.character(group)) {
    group <- read.csv(group, header = FALSE, row.names = NULL)
    colnames(group) <- c("Sample", "Group")
} else if (is.list(group)) {
    group <- do_call(
        rbind,
        lapply(names(group), function(n) data.frame(Sample = group[[n]], Group = n))
    )
} else if (!is.null(group)) {
    stop(paste0("Invalid group: ", paste0(group, collapse = ", ")))
}

logger$info("Reading and merging metrics for each sample ...")
metrics <- NULL
for (indir in indirs) {
    for (sens in sensitivity) {
        sample <- basename(indir)
        logger$debug("- Reading metrics for sample: ", sample)

        metric <- read.csv(
            file.path(indir, "metrics", paste0("sensitivity_", sens), "metrics_summary.csv"),
            header = FALSE, row.names = 1, check.names = FALSE)
        colnames(metric) <- c("Value")
        metric <- as.data.frame(t(metric))
        metric$Sample <- sample
        metric$Sensitivity <- sens

        if (!is.null(metrics)) {
            missing_cols <- setdiff(colnames(metrics), colnames(metric))
            if (length(missing_cols) > 0) {
                logger$warn('- Missing columns: {paste0(missing_cols, collapse = ", ")} in sample: {sample} (sensitivity={sens})')
                metric[missing_cols] <- NA
            }
            missing_cols <- setdiff(colnames(metric), colnames(metrics))
            if (length(missing_cols) > 0) {
                logger$warn('- Missing columns: {paste0(missing_cols, collapse = ", ")} in samples before {sample} (sensitivity={sens})')
                metrics[missing_cols] <- NA
            }
        }
        metrics <- rbind(metrics, metric)
    }
}

if (is.null(metrics)) {
    stop("No samples found, check the input directories.")
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
    list(kind = "descr", content = glue::glue("Metrics for all samples (sensitivity = {paste0(sensitivity, collapse = ", ")})")),
    list(kind = "table", src = file.path(outdir, "metrics.txt")),
    h1 = "Metrics of all samples"
)

METRIC_DESCR <- list(
    num_input_reads = "Total number of raw sequencing reads present in the input FASTQ files before any filtering.",
    num_reads_above_min_length = "Number of sequencing reads that pass the minimum length threshold during preprocessing. Reads shorter than the configured minimum length are discarded.",
    num_barcoded_reads = "Number of reads containing a valid cell barcode sequence that matches the expected barcode structure or whitelist.",
    num_mapped_reads = "Number of reads that successfully aligned to the reference transcriptome or genome after barcode processing.",
    num_mapped_reads_in_cells = "Number of mapped reads that are assigned to barcodes identified as true cells after cell-calling (i.e., excluding background or empty droplets).",
    mapping_pct_txome = "Percentage of reads that mapped to annotated transcriptomic regions (e.g., exonic regions consistent with gene models).",
    mapping_pct_genome = "Percentage of reads that mapped anywhere to the reference genome, including exonic, intronic, and intergenic regions.",
    num_cell_barcodes = "Number of distinct cell barcodes identified as real cells based on barcode rank distribution and cell-calling criteria.",
    mean_reads_per_cell = "Average number of reads assigned to each identified cell barcode. Typically calculated as total reads assigned to cells divided by the number of cells.",
    reads_in_cells_pct = "Percentage of mapped reads that are associated with identified cell barcodes rather than background barcodes.",
    duplication_rate = "Fraction of reads that are PCR duplicates, defined as reads sharing the same combination of cell barcode, UMI, and gene assignment. Reflects library complexity.",
    sequencing_saturation = "Fraction of reads originating from already-observed unique molecules (UMIs). Calculated as the proportion of duplicate reads among confidently mapped reads with valid barcodes and UMIs.",
    total_transcripts_in_cells = "Total number of deduplicated transcript molecules (UMIs) detected across all identified cells.",
    median_transcripts_in_cells = "Median number of deduplicated transcript molecules (UMIs) detected per identified cell.",
    num_genes_expressed_in_cells = "Total number of distinct genes with at least one detected transcript (UMI) across all identified cells.",
    median_genes_in_cells = "Median number of genes with at least one detected transcript (UMI) per identified cell.",
    pct_mito_in_cells = "Percentage of transcript counts in identified cells that are derived from mitochondrial genes. Commonly used as a cell quality control metric."
)

logger$info("Plotting metrics ...")
for (metric in colnames(metrics)) {
    if (metric == "Sample" || metric == "Sensitivity") { next }
    logger$info("- {metric}")

    reporter$add(
        list(
            kind = "descr",
            content = METRIC_DESCR[[metric]] %||% paste0("Metric: ", metric)
        ),
        h1 = metric
    )

    # barplot
    if (length(sensitivity) > 1) {
        p <- BarPlot(metrics, x = "Sample", y = metric, x_text_angle = 90, facet_by = "Sensitivity")
    } else {
        p <- BarPlot(metrics, x = "Sample", y = metric, x_text_angle = 90)
    }
    figfile <- file.path(outdir, paste0(slugify(metric), ".barplot"))
    save_plot(p, prefix = figfile)

    reporter$add(
        list(src = paste0(figfile, ".png"), name = "By Sample"),
        ui = "table_of_images",
        h1 = metric
    )

    if (is.null(group)) { next }
    # boxplot, if group is provided
    # group: Sample, Group
    pdata <- group %>%
        left_join(metrics, by = "Sample") %>%
        mutate(Group = factor(Group, levels = unique(Group)))

    if (length(sensitivity) > 1) {
        p <- BoxPlot(pdata, x = "Group", y = metric, x_text_angle = 90, facet_by = "Sensitivity")
    } else {
        p <- BoxPlot(pdata, x = "Group", y = metric, x_text_angle = 90)
    }
    figfile <- file.path(outdir, paste0(slugify(metric), ".boxplot"))
    save_plot(p, prefix = figfile)

    reporter$add(
        list(src = paste0(figfile, ".png"), name = "By Group"),
        ui = "table_of_images",
        h1 = metric
    )
}

reporter$save(joboutdir)
