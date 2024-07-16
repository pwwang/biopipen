
{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}

library(rlang)
library(ggplot2)
library(plotROC)

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}
joboutdir <- {{job.outdir | r}}
noids <- {{envs.noids | r}}
pos_label <- {{envs.pos_label | r}}
ci <- {{envs.ci | r}}
devpars <- {{envs.devpars | r}}
show_auc <- {{envs.show_auc | r}}
args <- {{envs.args | r: todot="-"}}
style_roc_args <- {{envs.style_roc | r: todot="-"}}
if (!is.null(style_roc_args$theme)) {
    style_roc_args$theme <- eval(parse(text=style_roc_args$theme))
}

data <- read.table(infile, header=TRUE, sep="\t", row.names = NULL, check.names = FALSE, stringsAsFactors=FALSE)
if (!noids) {
    data <- data[, -1]
}

# Normalize the first column (labels) into 0 and 1.
# If they are not 0/1, use pos_label to determine the positive class.
label_col <- colnames(data)[1]
if (is.character(data[[label_col]])) {
    data[[label_col]] <- as.numeric(data[[label_col]] == pos_label)
}

models <- colnames(data)[2:ncol(data)]

if (length(models) > 1) {
    # pivot longer the models, and put the model names into the column 'model'
    data <- melt_roc(data, label_col, colnames(data)[2:ncol(data)])
} else {
    data <- data.frame(
        D = data[[label_col]],
        M = data[[models]],
        name = rep(models, nrow(data))
    )
}

# Plot the ROC curve
p <- ggplot(data, aes(d = D, m = M, color = name))

if (isTRUE(ci)) {
    p <- p + do.call(geom_rocci, args)
} else {
    p <- p + do.call(geom_roc, args)
}

p <- p + do.call(style_roc, style_roc_args)
p <- p + scale_color_biopipen()

if (length(models) > 1) {
    p <- p + theme(legend.title = element_blank())
} else {
    p <- p + theme(legend.position = "none")
}

aucs = calc_auc(p)
write.table(aucs, file=file.path(joboutdir, "aucs.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

if (show_auc) {
    aucs = split(aucs$AUC, aucs$name)
    if (length(aucs) > 1) {
        # Add AUC values to the legend items
        p <- p +
            scale_color_manual(
                values = pal_biopipen()(length(models)),
                labels = sapply(models, function(m) paste(m, " (AUC =", round(aucs[[m]], 2), ")")),
                breaks = models)
    } else {
        p <- p +
            geom_text(
                x = 0.8, y = 0.2, label = paste("AUC =", round(unlist(aucs), 2)),
                color = "black", size = 4)
    }
}

devpars$filename <- outfile
do.call(png, devpars)
print(p)
dev.off()
