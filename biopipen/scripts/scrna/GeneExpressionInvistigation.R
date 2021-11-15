library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggprism)
library(ComplexHeatmap)

srtobjfile = {{in.srtobj | quote}}
genefile = {{in.genefile | quote}}
hmgenefile = {{in.hmgenefile | quote}}
groupfile = {{in.groupfile | quote}}
outdir = {{out.outdir | quote}}
config = {{in.configfile | read | toml_loads | r}}

if (length(config) == 0) {
    stop("No configuration specified.")
}

sobj = readRDS(srtobjfile)
genes = read.table(genefile, header=F, row.names=NULL, sep="\t", check.names=F)
hmgenes = read.table(hmgenefile, header=F, row.names=NULL, sep="\t", check.names=F)
if (ncol(genes) == 1) {
    colnames(genes) = c("Gene")
    genes$Name = genes$Gene
} else {
    colnames(genes) = c("Gene", "Name")
}
if (ncol(hmgenes) == 1) {
    colnames(hmgenes) = c("Gene")
    hmgenes$Name = hmgenes$Gene
} else {
    colnames(hmgenes) = c("Gene", "Name")
}

groups = read.table(groupfile, row.names=NULL, header=T, sep="\t", check.names = F)
groups = groups %>% rowwise() %>%
    mutate(across(-1, ~ strsplit(.x, ";", fixed=TRUE))) %>%
    as.data.frame()


allgenes = bind_rows(genes, hmgenes) %>% distinct(Gene, .keep_all = TRUE)

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
)[allgenes$Gene,]
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
    if (plottype == "boxplot") {
        cols = plotpms$ncol
        if (is.null(cols)) {
            cols = 3
        } else {
            plotpms$ncol = NULL
        }
        boxplotdata = plotdata %>%
            filter(Gene %in% genes$Name) %>%
            mutate(Gene = factor(Gene, levels=as.character(genes$Name)))

        p = ggplot(boxplotdata) +
            geom_boxplot(aes(x=Group, y=Log_Expression, fill=Group)) + facet_wrap(~Gene, ncol=cols) +
            theme_prism(axis_text_angle = 90) + theme(legend.position = "none") +
            xlab("")

        pngfile = file.path(outdir, paste0("boxplot.png"))
        plotpms$filename = pngfile
        do.call(png, plotpms)
        print(p)
        dev.off()

    } else {
        hmplotdata = plotdata %>%
            filter(Gene %in% hmgenes$Name) %>%
            pivot_wider(
                Group,
                names_from = Gene,
                values_from = Log_Expression,
                values_fn = mean
            ) %>%
            select(Group, all_of(hmgenes$Name)) %>%
            column_to_rownames("Group")

        p = Heatmap(
            hmplotdata,
            name = "Log_Expression",
            row_names_side = "left",
            # row_dend_side = "right",
            row_names_gp = gpar(fontsize = 12),
            row_names_max_width = max_text_width(
                rownames(hmplotdata),
                gp = gpar(fontsize = 12)
            ),
            # row_dend_width = unit(2, "cm"),
            rect_gp = gpar(col = "#DDDDDD", lwd = 1),
            cluster_columns = FALSE,
            cluster_rows = FALSE
        )
        pngfile = file.path(outdir, paste0("heatmap.png"))
        plotpms$filename = pngfile
        do.call(png, plotpms)
        print(p)
        dev.off()
    }
}

