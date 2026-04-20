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
