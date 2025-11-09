library(rlang)
library(Seurat)
library(scplotter)
library(biopipen.utils)

adfile <- {{in.adfile | r}}
outfile <- {{out.outfile | r}}
dotplot_check <- {{envs.dotplot_check | r}}
outdir <- dirname(outfile)
assay <- {{envs.assay | r}}
ident <- {{envs.ident | r}}

log <- get_logger()

ConvertAnnDataToSeurat(adfile, outfile = outfile, assay = assay, ident = ident, log = log)

if (!isFALSE(dotplot_check)) {
    log$info("Reading Seurat object ...")
    sobj <- read_obj(outfile)

    log$info("Checking dotplot ...")
    dotfig <- file.path(outdir, "dotplot.png")
    if (isTRUE(dotplot_check)) {
        vobj <- FindVariableFeatures(
            sobj, selection.method = "vst", nfeatures = 2000)
        dotplot_check <- head(VariableFeatures(vobj), 10)
    } else if (is.character(dotplot_check)) {
        dotplot_check <- trimws(strsplit(dotplot_check, ",")[[1]])
    }
    p <- FeatureStatPlot(
        sobj, features = dotplot_check, plot_type = "dot",
        assay = assay
    )
    res = 70
    height <- attr(p, "height") * res
    width <- attr(p, "width") * res
    png(dotfig, width = width, height = height, res = res)
    print(p)
    dev.off()
}
