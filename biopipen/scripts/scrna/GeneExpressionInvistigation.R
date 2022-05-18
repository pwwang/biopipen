source("{{biopipen_dir}}/utils/io.R")
source("{{biopipen_dir}}/utils/plot.R")
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(rlang)
library(ggplot2)
library(ggprism)
library(ComplexHeatmap)

srtobjfile = {{in.srtobj | quote}}
genefile = {{in.genefile | r}}
outdir = {{out.outdir | quote}}
gopts = {{envs.gopts | r}}
{% if in.configfile %}
config = {{in.configfile | toml_load | r}}
{% set config = in.configfile | toml_load %}
{% else %}
config = {{envs.config | r}}
{% set config = envs.config %}
{% endif %}

sobj = readRDS(srtobjfile)
genes = read.table.opts(genefile, gopts)
if (ncol(genes) == 1) {
    genes$.Name = genes[[1]]
}
colnames(genes) = c("Gene", "Name")

if (!is.null(config$mutaters)) {
    expressions = list()
    for (key in names(config$mutaters)) {
        expressions[[key]] = parse_expr(config$mutaters[[key]])
    }
    sobj@meta.data = mutate(sobj@meta.data, !!!expressions)
}

if (!is.null(config$subset)) {
    sobj = subset(sobj, subset = {{config.subset}})
}

DefaultAssay(sobj) <- "RNA"
sobj = NormalizeData(sobj)

exprs = as.data.frame(
    GetAssayData(sobj, slot = "data", assay = "RNA")
)[genes$Gene,,drop=F]
rownames(exprs) = genes$Name
exprs = rownames_to_column(exprs, "Gene")

plot_heatmap = function(plotconf, outfile) {
    plotdata = exprs |>
        pivot_longer(
            names(exprs)[2:ncol(exprs)],
            names_to = "Barcode",
            values_to = "Log_Expression"
        )
    metadata = sobj@meta.data[plotdata$Barcode,,drop=F]
    plotdata = cbind(plotdata, metadata)
    plotdata = plotdata |>
        group_by(Gene, !!sym(config$groupby)) |>
        summarise(Log_Expression = mean(Log_Expression)) |>
        pivot_wider(names_from = config$groupby, values_from = "Log_Expression") |>
        column_to_rownames("Gene")

    given_genes = rownames(plotdata)
    plotdata = plotdata[complete.cases(plotdata),,drop=FALSE]
    invalid_genes = setdiff(given_genes, rownames(plotdata))
    if (length(invalid_genes) > 0) {
        warning(
            paste(
                "The following genes were not found in the data:",
                invalid_genes
            )
        )
    }

    devpars = list(res=plotconf$res, width=plotconf$width, height=plotconf$height)
    plotconf$res = NULL
    plotconf$width = NULL
    plotconf$height = NULL

    for (name in names(plotconf)) {
        plotconf[[name]] = parse_expr(plotconf[[name]])
    }

    plotHeatmap(
        plotdata,
        plotconf,
        devpars = devpars,
        outfile = outfile
    )
}

plot_boxplot = function(plotconf, outfile) {
    plotdata = exprs |>
        pivot_longer(
            names(exprs)[2:ncol(exprs)],
            names_to = "Barcode",
            values_to = "Log_Expression"
        )
    metadata = sobj@meta.data[plotdata$Barcode,,drop=F]
    plotdata = cbind(plotdata, metadata)

    cols = if (is.null(plotconf$ncol)) 3 else plotconf$ncol
    p = ggplot(plotdata) +
        geom_boxplot(aes_string(x=config$groupby, y="Log_Expression", fill=config$groupby)) +
        facet_wrap(~Gene, ncol=cols) +
        theme_prism(axis_text_angle = 90) + theme(legend.position = "none") +
        xlab("")

    devpars = list(filename = outfile, res = plotconf$res, width = plotconf$width, height = plotconf$height)
    do.call(png, devpars)
    print(p)
    dev.off()
}


for (plottype in names(config$plots)) {
    plotconf = config$plots[[plottype]]
    if (plottype == "heatmap") {
        plotfile = file.path(outdir, "heatmap.png")
        plot_heatmap(plotconf, plotfile)
    } else if (plottype == "boxplot") {
        plotfile = file.path(outdir, "boxplot.png")
        plot_boxplot(plotconf, plotfile)
    } else {
        stop(paste("Unknown plot type:", plottype))
    }
}
