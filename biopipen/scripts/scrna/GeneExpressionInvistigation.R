source("{{biopipen_dir}}/utils/io.R")
source("{{biopipen_dir}}/utils/plot.R")
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggprism)
library(ComplexHeatmap)

{% addfilter compile_config %}
def compile_config(tomldict):
    out = {}
    for key, val in tomldict.items():
        if isinstance(val, str) and val[:2].lower() == "r:":
            # modify "r:" to "@r:" to preserve the expression from "r()"
            # since we don't have the context variable ready here
            out[key] = f"@{val}"
        elif isinstance(val, dict):
            out[key] = compile_config(val)
        else:
            out[key] = val
    return out
{% endaddfilter %}

srtobjfile = {{in.srtobj | quote}}
groupfile = {{in.groupfile | quote}}
genefiles = {{in.genefiles | r}}
outdir = {{out.outdir | quote}}
envs = {{envs | r}}
config = {{in.configfile | read | toml_loads | compile_config | r}}
for (name in names(envs)) {
    if (is.null(config[[name]])) {
        config[[name]] = envs[[name]]
    }
}
if (is.null(config$name)) {
    config$name = {{in.groupfile | stem | r}}
}

sobj = readRDS(srtobjfile)
genes = lapply(genefiles, function(x) {
    out = read.table.opts(x, config$gopts)
    if (ncol(out) == 1) {
        out$.Name = out[[1]]
    }
    colnames(out) = c("Gene", "Name")
    return(out)
})
allgenes = do.call(bind_rows, genes) %>% distinct(Gene, .keep_all = T)

groups = read.table(groupfile, row.names=NULL, header=T, sep="\t", check.names = F)
groups = groups %>% rowwise() %>%
    mutate(across(-1, ~ strsplit(.x, ";", fixed=TRUE))) %>%
    as.data.frame()

samples = config$target
cells = c()
# cells = paste(samples[1], unlist(groups[[ samples[1] ]]), sep="_")
# merge seurat object and add cell ids

# y = c()
if (length(samples) > 1) {
    for (i in 1:length(samples)) {
        # y = c(y, seurat_obj[[samples[i]]])
        cells = c(
            cells,
            paste(samples[i], unlist(groups[[ samples[i] ]]), sep="_")
        )
    }
}
# if (length(y) == 0) {
#     sobj = seurat_obj[[samples[1]]]
#     sobj = RenameCells(sobj, add.cell.id = samples[1])
# } else {
#     sobj = merge(seurat_obj[[samples[1]]], y, add.cell.ids = samples)
# }
sobj = subset(sobj, cells = cells)
# already Normalized by SeuratPreparing
# sobj = NormalizeData(sobj, normalization.method = "LogNormalize")

DefaultAssay(sobj) <- "RNA"
sobj = NormalizeData(sobj)

exprs = as.data.frame(
    GetAssayData(sobj, slot = "data", assay = "RNA")
)[allgenes$Gene,,drop=F]
rownames(exprs) = allgenes$Name

exprdata = list()
i = 1
for (group in unique(groups[[1]])) {
    barcodes = groups %>% filter(.[[1]] == group)
    bcodes = c()
    for (sample in samples) {
        bcodes = c(
            bcodes,
            paste(sample, unlist(barcodes[[sample]]), sep="_")
        )
    }
    bcodes = intersect(bcodes, colnames(exprs))

    if (length(bcodes) == 0) {
        exprdata[[i]] = NULL
    } else {
        exprdata[[i]] = exprs[, bcodes, drop=F] %>%
            as.data.frame() %>%
            rownames_to_column("Gene") %>%
            pivot_longer(-"Gene", names_to="Barcode", values_to="Log_Expression") %>%
            mutate(Group=group)
    }
    i = i + 1
}
plotdata = do.call(bind_rows, exprdata)

for (plottype in names(config$plots)) {
    plotpms = config$plots[[plottype]]
    if (is.null(plotpms)) { plotpms = list() }
    if (is.null(plotpms$use)) {
        plotgenes = genes[[1]]
    } else {
        plotgenes = genes[[plotpms$use]]
        plotpms$use = NULL
    }
    myplotdata = plotdata %>%
        filter(Gene %in% plotgenes$Name) %>%
        mutate(Gene = factor(Gene, levels=as.character(plotgenes$Name)))

    if (plottype == "boxplot") {
        cols = plotpms$ncol
        if (is.null(cols)) {
            cols = 3
        } else {
            plotpms$ncol = NULL
        }
        p = ggplot(myplotdata) +
            geom_boxplot(aes(x=Group, y=Log_Expression, fill=Group)) +
            facet_wrap(~Gene, ncol=cols) +
            theme_prism(axis_text_angle = 90) + theme(legend.position = "none") +
            xlab("")

        pngfile = file.path(outdir, paste0("boxplot.png"))
        plotpms$filename = pngfile
        do.call(png, plotpms)
        print(p)
        dev.off()

    } else {
        hmdata = myplotdata %>%
            pivot_wider(
                Group,
                names_from = Gene,
                values_from = Log_Expression,
                values_fn = mean
            ) %>%
            select(Group, all_of(plotgenes$Name)) %>%
            column_to_rownames("Group")

        devpars = list()
        devpars$res = plotpms$res
        devpars$height = plotpms$height
        devpars$width = plotpms$width
        plotpms$res = NULL
        plotpms$height = NULL
        plotpms$width = NULL

        for (name in names(plotpms)) {
            if (is.character(plotpms[[name]]) &&
                startsWith(plotpms[[name]], "@r:")) {
                plotpms[[name]] = eval(
                    parse(text=substring(plotpms[[name]], 4))
                )
            }
        }

        pngfile = file.path(outdir, paste0("heatmap.png"))
        plotHeatmap(
            hmdata,
            plotpms,
            devpars = devpars,
            outfile = pngfile
        )

    }
}

