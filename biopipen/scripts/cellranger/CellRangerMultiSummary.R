library(rlang)
library(dplyr)
library(tidyr)
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

logger$info("Reading and merging metrics for each run/sample ...")
metrics <- NULL
for (indir in indirs) {
    run_name <- basename(indir)
    per_sample_dir <- file.path(indir, "outs", "per_sample_outs")
    if (!dir.exists(per_sample_dir)) {
        logger$warn("- per_sample_outs not found for run: ", run_name, ", skipping.")
        next
    }
    sample_ids <- list.dirs(per_sample_dir, full.names = FALSE, recursive = FALSE)
    for (sample_id in sample_ids) {
        sample_name <- paste0(run_name, "/", sample_id)
        logger$debug("- Reading metrics for: ", sample_name)
        metrics_file <- file.path(per_sample_dir, sample_id, "metrics_summary.csv")
        if (!file.exists(metrics_file)) {
            logger$warn("  - metrics_summary.csv not found, skipping.")
            next
        }
        metric <- read.csv(metrics_file, header = TRUE, row.names = NULL, check.names = FALSE)
        metric$Sample <- sample_name

        if (!is.null(metrics)) {
            missing_cols <- setdiff(colnames(metrics), colnames(metric))
            if (length(missing_cols) > 0) {
                logger$warn('  - Missing columns: {paste0(missing_cols, collapse = ", ")} in: {sample_name}')
                metric[missing_cols] <- NA
            }
            missing_cols <- setdiff(colnames(metric), colnames(metrics))
            if (length(missing_cols) > 0) {
                logger$warn('  - Missing columns: {paste0(missing_cols, collapse = ", ")} in earlier samples')
                metrics[missing_cols] <- NA
            }
        }
        metrics <- rbind(metrics, metric)
    }
}

if (is.null(metrics)) {
    stop("No samples found, check the input directories.")
}

# cellranger "multi" writes metrics_summary.csv in long format (one row per metric) and each per-sample
# file also repeats the run/library-level rows (Grouped By = "Physical library ID", "Fastq ID",
# "Probe barcode ID", ...). Only the Category == "Cells" rows describe the sample itself. Pivot them into
# the one-row-per-sample / one-column-per-metric layout the rest of this script expects (the layout a
# CellRangerCount metrics_summary.csv already has).
metrics <- metrics %>%
    filter(Category == "Cells") %>%
    dplyr::select(Sample, MetricName = `Metric Name`, MetricValue = `Metric Value`) %>%
    pivot_wider(names_from = MetricName, values_from = MetricValue, values_fn = first)

# Percent columns come in two forms: "40.8%" (count pipelines) and plain fractions such as 0.408
# (cellranger 10 "multi" metrics_summary.csv). Detect whichever form this file uses; a file with any
# "%" value keeps the original handling exactly, so count summaries are unaffected.
pct_symbol <- function(v) is.character(v) && any(grepl("%", v[!is.na(v)]))
pct_fraction <- function(v) {
    if (!is.character(v)) return(FALSE)
    v <- v[!is.na(v) & v != ""]
    if (length(v) == 0 || any(grepl("%", v))) return(FALSE)
    nums <- suppressWarnings(as.numeric(v))
    !any(is.na(nums)) && all(nums >= 0 & nums <= 1)
}
has_symbols <- any(sapply(metrics, pct_symbol))
detect <- if (has_symbols) pct_symbol else pct_fraction
percent_columns <- setdiff(colnames(metrics)[sapply(metrics, detect)], "Sample")
percent_symbol <- if (has_symbols) percent_columns else character(0)
percent_fraction <- if (has_symbols) character(0) else percent_columns
# Remove % / scale fractions, then label those columns as (%). Guard the empty case: rename_with()
# aborts when the selection is empty (metrics with no percent-like column at all).
if (length(percent_columns) > 0) {
    metrics <- metrics %>%
        mutate(across(all_of(percent_symbol), ~as.numeric(gsub("%", "", .x)))) %>%
        mutate(across(all_of(percent_fraction), ~as.numeric(.x) * 100)) %>%
        rename_with(.fn = function(x) { paste0(x, " (%)") }, .cols = all_of(percent_columns))
}
metrics <- metrics %>%
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
    list(kind = "descr", content = "Metrics for all samples across all runs"),
    list(kind = "table", src = file.path(outdir, "metrics.txt")),
    h1 = "Metrics of all samples"
)

logger$info("Plotting metrics ...")
for (metric in colnames(metrics)) {
    if (metric == "Sample") { next }
    metric_name <- sub(" \\(%\\)$", "", metric)
    logger$info("- {metric_name}")

    reporter$add(
        list(
            kind = "descr",
            content = paste0("Metric: ", metric)
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
