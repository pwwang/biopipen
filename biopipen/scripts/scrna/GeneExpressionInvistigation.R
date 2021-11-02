library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggprism)

srtobjfile = {{in.srtobj | quote}}
genefile = {{in.genefile | quote}}
groupfile = {{in.groupfile | quote}}
outdir = {{out.outdir | quote}}
cases = {{envs.cases | r}}

if (length(cases) == 0) {
    stop("No `envs.cases` specified.")
}

seurat_obj = readRDS(srtobjfile)
genes = read.table(genefile, header=F, row.names=NULL, sep="\t", check.names=F)
if (ncol(genes) == 1) {
    colnames(genes) = c("Gene")
} else {
    colnames(genes) = c("Gene", "Name")
}

groups = read.table(groupfile, row.names=1, header=T, sep="\t", check.names = F)
gnames = rownames(groups)
groups = groups %>% rowwise() %>%
    mutate(across(everything(), ~ strsplit(.x, ";", fixed=TRUE))) %>%
    as.data.frame()
group_ns = groups
group_ns[is.na(group_ns)] = list(NULL)
group_ns = group_ns %>% mutate(across(everything(), lengths))
rownames(groups) = gnames
rownames(group_ns) = gnames

x = 0
for (case in names(cases)) {
    gdata_ns = group_ns %>% filter(eval(parse(text=cases[[case]]$condition)))
    gdata = groups[rownames(gdata_ns),,drop=F]
    samples = cases[[case]]$target
    cells = paste(samples[1], unlist(gdata[[ samples[1] ]]), sep="_")
    # merge seurat object and add cell ids

    y = c()
    if (length(samples) > 1) {
        for (i in 2:length(samples)) {
            y = c(y, seurat_obj[[samples[i]]])
            cells = c(
                cells,
                paste(samples[i], unlist(gdata[[ samples[i] ]]), sep="_")
            )
        }
    }
    if (length(y) == 0) {
        sobj = seurat_obj[[samples[1]]]
        sobj = RenameCells(sobj, add.cell.id = samples[1])
    } else {
        sobj = merge(seurat_obj[[samples[1]]], y, add.cell.ids = samples)
    }
    sobj = subset(sobj, cells = cells)
    sobj = NormalizeData(sobj, normalization.method = "LogNormalize")

    exprs = as.data.frame(
        GetAssayData(sobj, slot = "data", assay = "RNA")
    )[as.character(genes$Gene),]
    if (ncol(genes) > 1) {
        rownames(exprs) = as.character(genes$Name)
    }

    exprdata = list()
    i = 1
    for (group in rownames(gdata)) {
        barcodes = gdata[group, samples, drop=FALSE]
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
    plotdata = do.call(bind_rows, exprdata) %>%
        mutate(Gene = factor(Gene, levels=as.character(genes$Name)))

    cols = cases[[case]]$ncol
    if (is.null(cols)) {
        cols = 3
    }
    p = ggplot(plotdata) +
        geom_boxplot(aes(x=Group, y=Log_Expression, fill=Group)) + facet_wrap(~Gene, ncol=cols) +
        theme_prism(axis_text_angle = 90) + theme(legend.position = "none") +
        xlab("") +
        ggtitle(case)

    slug = cases[[case]]$slug
    if (is.null(slug)) {
        slug = x
    }
    cat(case, file = file.path(outdir, paste0(slug, ".title")))
    pngfile = file.path(outdir, paste0(slug, ".png"))
    if (is.null(cases[[case]]$devpars)) {
        cases[[case]]$devpars = list()
    }
    pngparams = c(pngfile, cases[[case]]$devpars)
    do.call(png, pngparams)
    print(p)
    dev.off()

    x = x + 1
}

