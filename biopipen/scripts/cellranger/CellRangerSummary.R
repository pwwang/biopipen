library(rlang)
library(dplyr)
library(plotthis)
library(biopipen.utils)

indirs <- {{in.indirs | r}}
outdir <- {{out.outdir | r}}
joboutdir <- {{job.outdir | r}}
group <- {{envs.group | r}}

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

cellranger_type <- NULL
logger$info("Reading and merging metrics for each sample ...")
metrics <- NULL
for (indir in indirs) {
    sample <- basename(indir)
    logger$debug("- Reading metrics for sample: ", sample)
    metric <- read.csv(
        file.path(indir, "outs", "metrics_summary.csv"),
        header = TRUE, row.names = NULL, check.names = FALSE)
    metric$Sample <- sample
    sample_cellranger_type <- ifelse(
        file.exists(file.path(indir, "outs", "clonotypes.csv")),
        "vdj",
        "count"  # support more types in the future
    )
    cellranger_type <- cellranger_type %||% sample_cellranger_type
    if (cellranger_type != sample_cellranger_type) {
        stop("Multiple types of CellRanger output detected. Should be either count or vdj.")
    }
    if (!is.null(metrics)) {
        missing_cols <- setdiff(colnames(metrics), colnames(metric))
        if (length(missing_cols) > 0) {
            logger$warn('- Missing columns: {paste0(missing_cols, collapse = ", ")} in sample: {sample}')
            metric[missing_cols] <- NA
        }
        missing_cols <- setdiff(colnames(metric), colnames(metrics))
        if (length(missing_cols) > 0) {
            logger$warn('- Missing columns: {paste0(missing_cols, collapse = ", ")} in samples before {sample}')
            metrics[missing_cols] <- NA
        }
    }
    metrics <- rbind(metrics, metric)
}

if (is.null(metrics)) {
    stop("No samples found, check the input directories.")
}

percent_columns <- sapply(colnames(metrics), function(x) {
    is.character(metrics[[x]]) && grepl("%", metrics[[x]][1]) && x != "Sample"
})
percent_columns <- colnames(metrics)[percent_columns]
# Remove %
metrics <- metrics %>%
    mutate(across(all_of(percent_columns), ~as.numeric(gsub("%", "", .x)))) %>%
    rename_with(.fn = function(x) { paste0(x, " (%)") }, .cols = percent_columns) %>%
    mutate(across(-Sample, ~as.numeric(gsub(",", "", .x))))

# Save metrics
write.table(
    metrics,
    file.path(outdir, "metrics.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

reporter$add(
    list(kind = "descr", content = "Metrics for all samples"),
    list(kind = "table", src = file.path(outdir, "metrics.txt")),
    h1 = "Metrics of all samples"
)

if (cellranger_type == "vdj") {
    METRIC_DESCR = list(
        `Estimated Number of Cells` = "The number of barcodes estimated to correspond to GEMs containing cells. See VDJ Cell Calling Algorithm.",
        `Mean Read Pairs per Cell` = "Number of input read pairs divided by the estimated number of cells.",
        `Number of Cells With Productive V-J Spanning Pair` = "Number of cell barcodes for which at least one productive sequence was found for each of TRA and TRB (or heavy and light chains, for Ig).",
        `Number of Read Pairs` = "Total number of read pairs that were assigned to this library in demultiplexing.",
        `Valid Barcodes` = "Fraction of reads with barcodes that match the whitelist after barcode correction.",
        `Q30 Bases in Barcode` = "Fraction of cell barcode bases with Q-score greater than or equal to 30.",
        `Q30 Bases in RNA Read 1` = "Fraction of read 1 bases with Q-score greater than or equal to 30. (Likewise for read 2.)",
        `Q30 Bases in Sample Index` = "Fraction of sample index bases with Q-score greater than or equal to 30.",
        `Q30 Bases in UMI` = "Fraction of UMI bases with Q-score ≥ 30.",
        `Reads Mapped to Any V(D)J Gene` = "Fraction of reads that partially or wholly map to a V(D)J gene segment.",
        `Reads Mapped to TRA` = "Fraction of reads that map partially or wholly to a TRA gene segment.",
        `Mean Used Read Pairs per Cell` = "Mean number of read pairs used in assembly per cell barcode. These reads must have a cell barcode, map to a V(D)J gene, and have a UMI with sufficient read support, counted after subsampling.",
        `Fraction Reads in Cells` = "Number of reads with cell barcodes divided by the number of reads with valid barcodes.",
        `Median TRA UMIs per Cell` = "Median number of UMIs assigned to a TRA contig per cell.",
        `Paired Clonotype Diversity` = "Effective diversity of the paired clonotypes, computed as the Inverse Simpson Index of the clonotype frequencies. A value of 1 indicates a minimally diverse sample - only one distinct clonotype was detected. A value equal to the estimated number of cells indicates a maximally diverse sample.",
        `Cells With TRA Contig` = "Fraction of cell barcodes with at least one TRA contig annotated as a full or partial V(D)J gene.",
        `Cells With CDR3-annotated TRA Contig` = "Fraction of cell barcodes with at least one TRA contig where a CDR3 was detected.",
        `Cells With V-J Spanning Contig` = "Fraction of cell barcodes with at least one full-length contig.",
        `Cells With V-J Spanning TRA Contig` = "Fraction of cell barcodes with at least one full-length TRA contig.",
        `Cells With Productive TRA Contig` = "Fraction of cell barcodes with at least one full-length TRA contig that is productive.",
        `Cells With Productive V-J Spanning Pair` = "Fraction of cell barcodes with at least one contig for each chain of the receptor pair that is productive."
    )
} else {
    METRIC_DESCR = list(
        `Estimated Number of Cells` = "The number of barcodes associated with cell-containing partitions, estimated from the barcode UMI count distribution.",
        `Mean Reads per Cell` = "The total number of reads divided by the estimated number of cells.",
        `Median Genes per Cell` = "Median number of read pairs sequenced from the cells assigned to this sample. In case of multiplexing, only cell-associated barcodes assigned exactly one CMO can be assigned to a sample.",
        `Number of Reads` = "Total number of sequencing reads.",
        `Valid Barcodes` = "Fraction of reads with cell-barcodes that match the whitelist.",
        `Sequencing Saturation` = 'Fraction of reads originating from an already-observed UMI. This is a function of library complexity and sequencing depth. More specifically, this is a ratio where: the denominator is the number of confidently-mapped reads with a valid cell-barcode and valid UMI, and the numerator is the subset of those reads that had a non-unique combination of (cell-barcode, UMI, gene). This metric was called "cDNA PCR Duplication" in versions of Cell Ranger prior to 1.2.',
        `Q30 Bases in Barcode` = "Fraction of bases with Q-score at least 30 in the cell barcode sequences. This is the i7 index (I1) read for the Single Cell 3' v1 chemistry and the R1 read for the Single Cell 3' v2 chemistry.",
        `Q30 Bases in RNA` = "Fraction of bases with Q-score at least 30 in the RNA read sequences. This is Illumina R1 for the Single Cell 3' v1 chemistry and Illumina R2 for the Single Cell 3' v2 chemistry.",
        `Q30 Bases in UMI` = "Fraction of bases with Q-score at least 30 in the UMI sequences. This is the R2 read for the Single Cell 3' v1 chemistry and the R1 read for the Single Cell 3' v2 chemistry.",
        `Reads Mapped to Genome` = "Fraction of reads that mapped to the genome.",
        `Reads Mapped Confidently to Genome` = "Fraction of reads that mapped uniquely to the genome. If a read mapped to exonic loci from a single gene and also to non-exonic loci, it is considered uniquely mapped to one of the exonic loci.",
        `Reads Mapped Confidently to Intergenic Regions` = "Fraction of reads that mapped to the intergenic regions of the genome with a high mapping quality score as reported by the aligner.",
        `Reads Mapped Confidently to Intronic Regions` = "Fraction of reads that mapped to the intronic regions of the genome with a high mapping quality score as reported by the aligner.",
        `Reads Mapped Confidently to Exonic Regions` = "Fraction of reads that mapped to the exonic regions of the genome with a high mapping quality score as reported by the aligner.",
        `Reads Mapped Confidently to Transcriptome` = "Fraction of reads that mapped to a unique gene in the transcriptome with a high mapping quality score as reported by the aligner. The read must be consistent with annotated splice junctions when include-introns=false. These reads are considered for UMI counting.",
        `Reads Confidently Mapped Antisense` = "Fraction of reads confidently mapped to the transcriptome, but on the opposite strand of their annotated gene. A read is counted as antisense if it has any alignments that are consistent with an exon of a transcript but antisense to it, and has no sense alignments.",
        `Total Genes Detected Median UMI Counts per Cell` = "The number of genes with at least one UMI count in any cell."
    )
}
logger$info("Plotting metrics ...")
for (metric in colnames(metrics)) {
    if (metric == "Sample") { next }
    metric_name <- sub(" \\(%\\)$", "", metric)
    logger$info("- {metric_name}")

    reporter$add(
        list(
            kind = "descr",
            content = METRIC_DESCR[[metric_name]] %||% paste0("Metric: ", metric)
        ),
        h1 = metric
    )

    # barplot
    p <- BarPlot(metrics, x = "Sample", y = metric, x_text_angle = 90)
    figfile <- file.path(outdir, paste0(slugify(metric), ".barplot.png"))
    png(figfile, height = 600, res = 100, width = max(nrow(metrics) * 30 + 200, 400))
    print(p)
    dev.off()

    reporter$add(
        list(src = figfile, name = "By Sample"),
        ui = "table_of_images",
        h1 = metric
    )

    if (is.null(group)) { next }
    # boxplot, if group is provided
    # group: Sample, Group
    pdata <- group %>%
        left_join(metrics, by = "Sample") %>%
        mutate(Group = factor(Group, levels = unique(Group)))

    p <- BoxPlot(pdata, x = "Group", y = metric, x_text_angle = 90)
    figfile <- file.path(outdir, paste0(slugify(metric), ".boxplot.png"))
    png(figfile, height = 600, res = 100, width = max(length(unique(pdata$Group)) * 30 + 200, 400))
    print(p)
    dev.off()

    reporter$add(
        list(src = figfile, name = "By Group"),
        ui = "table_of_images",
        h1 = metric
    )
}

reporter$save(joboutdir)
