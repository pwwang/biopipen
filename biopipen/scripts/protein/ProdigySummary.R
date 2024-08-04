{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}

library(rlang)
library(dplyr)
library(ggplot2)
library(ggprism)

theme_set(theme_prism())

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

log_info("Reading and merging metrics for each sample ...")
metrics <- NULL

for (infile in infiles) {
    sample <- sub("_prodigy$", "", basename(dirname(infile)))
    log_debug("- Reading metrics from {sample}")
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

add_report(
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
    log_info("Merging group information ...")
    metrics <- group %>%
        left_join(metrics, by = "Sample") %>%
        mutate(Group = factor(Group, levels = unique(Group)))
}

log_info("Plotting Prodigy metrics ...")
for (metric in names(METRIC_DESCR)) {
    log_info("- {metric}: {METRIC_DESCR[[metric]]}")

    add_report(
        list(
            kind = "descr",
            content = METRIC_DESCR[[metric]] %||% paste0("Metric: ", metric)
        ),
        h1 = metric
    )

    # barplot
    p <- ggplot(metrics, aes(x = Sample, y = !!sym(metric))) +
        geom_bar(stat = "identity", fill = "steelblue") +
        labs(x = "Sample", y = metric) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    figfile <- file.path(outdir, paste0(slugify(metric), ".barplot.png"))
    png(figfile, height = 600, res = 100, width = nrow(metrics) * 30 + 200)
    print(p)
    dev.off()

    add_report(
        list(src = figfile, name = "By Sample"),
        ui = "table_of_images",
        h1 = metric
    )

    if (is.null(group)) { next }
    # group: Sample, Group
    p <- ggplot(metrics, aes(x = Group, y = !!sym(metric))) +
        geom_boxplot(fill = "steelblue") +
        labs(x = "Group", y = metric) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    figfile <- file.path(outdir, paste0(slugify(metric), ".boxplot.png"))
    png(figfile, height = 600, res = 100, width = length(unique(metrics$Group)) * 30 + 200)
    print(p)
    dev.off()

    add_report(
        list(src = figfile, name = "By Group"),
        ui = "table_of_images",
        h1 = metric
    )
}

save_report(joboutdir)
