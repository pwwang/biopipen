library(rlang)
library(dplyr)
library(biopipen.utils)
library(plotthis)

infiles <- {{in.infiles | r}}
outdir <- {{out.outdir | r}}
joboutdir <- {{job.outdir | r}}
group <- {{envs.group | r}}

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

log <- get_logger()
reporter <- get_reporter()

log$info("Reading and merging metrics for each sample ...")
metrics <- NULL

for (infile in infiles) {
    sample <- sub("_prodigy$", "", basename(dirname(infile)))
    log$debug("- Reading metrics from {sample}")
    metric <- read.table(
        infile,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        row.names = NULL)
    metric$Sample <- sample
    metric <- metric %>% select(Sample, everything())
    if (is.null(metrics)) {
        metrics <- metric
    } else {
        metrics <- rbind(metrics, metric)
    }
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
    list(kind = "descr", content = "Metrics for all samples"),
    list(kind = "table", src = file.path(outdir, "metrics.txt")),
    h1 = "Metrics of all samples"
)

METRIC_DESCR = list(
    nIC = "No. of intermolecular contacts",
    nCCC = "No. of charged-charged contacts",
    nCPC = "No. of charged-polar contacts",
    nCAPC = "No. of charged-apolar contacts",
    nPPC = "No. of polar-polar contacts",
    nAPPC = "No. of apolar-polar contacts",
    nAPAPC = "No. of apolar-apolar contacts",
    pANISR = "Percentage of apolar NIS residues",
    pCNISR = "Percentage of charged NIS residues",
    BindingAffinity = "Predicted binding affinity (kcal.mol^-1)",
    DissociationConstant = "Predicted dissociation constant (M)"
)

if (!is.null(group)) {
    log$info("Merging group information ...")
    metrics <- group %>%
        left_join(metrics, by = "Sample") %>%
        mutate(Group = factor(Group, levels = unique(Group)))
}

log$info("Plotting Prodigy metrics ...")
for (metric in names(METRIC_DESCR)) {
    log$info("- {metric}: {METRIC_DESCR[[metric]]}")

    reporter$add(
        list(
            kind = "descr",
            content = METRIC_DESCR[[metric]] %||% paste0("Metric: ", metric)
        ),
        h1 = metric
    )

    p <- plotthis::BarPlot(
        x = "Sample",
        y = metric,
        x_text_angle = 90,
        fill = "Group",
        data = metrics
    )

    figfile <- file.path(outdir, paste0(slugify(metric), ".barplot.png"))
    height <- attr(p, "height") %||% 6
    width <- attr(p, "width") %||% (nrow(metrics) * .3 + 2)
    png(figfile, height = height * 100, res = 100, width = width * 100)
    print(p)
    dev.off()

    reporter$add(
        list(src = figfile, name = "By Sample"),
        ui = "table_of_images",
        h1 = metric
    )

    if (is.null(group)) { next }
    # group: Sample, Group
    p <- plotthis::BarPlot(
        data = metrics,
        x = "Group",
        y = metric,
        x_text_angle = 90
    )

    figfile <- file.path(outdir, paste0(slugify(metric), ".boxplot.png"))
    height <- attr(p, "height") %||% 6
    width <- attr(p, "width") %||% (length(unique(metrics$Group)) * 0.3 + 2)
    png(figfile, height = height * 100, res = 100, width = width * 100)
    print(p)
    dev.off()

    reporter$add(
        list(src = figfile, name = "By Group"),
        ui = "table_of_images",
        h1 = metric
    )
}

reporter$save(joboutdir)
