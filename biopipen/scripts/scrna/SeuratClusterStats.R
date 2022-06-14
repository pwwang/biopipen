source("{{biopipen_dir}}/utils/plot.R")
library(Seurat)
library(dplyr)
library(ggprism)

srtfile = {{in.srtobj | r}}
outdir = {{out.outdir | r}}
envs = {{envs | r}}

srtobj = readRDS(srtfile)


do_stats_ncells = function(df_ncells, odir, pms) {
    outfile = file.path(odir, "ncells.png")
    plotGG(
        df_ncells,
        "col",
        list(
            mapping = aes(
                x=seurat_clusters,
                y=nCellsPerCluster,
                fill=seurat_clusters
            )
        ),
        c(
            'ggtitle("Total # cells for each cluster")',
            'theme_prism()',
            'labs(x="Seurat Cluster", y="# Cells Per Cluster")'
        ),
        pms$devpars,
        outfile
    )
}


do_stats_ncellspersample = function(df_ncells, odir, pms) {
    outfile = file.path(odir, "ncellspersample.png")
    plotGG(
        df_ncells,
        "col",
        list(
            mapping = aes(
                x=seurat_clusters,
                y=nCellsPerCluster,
                fill=Sample
            ),
            position="stack"
        ),
        c(
            'ggtitle("Breakdown # cells from each samples for each cluster")',
            'theme_prism()',
            'labs(x="Seurat Cluster", y="# Cells Per Cluster From Samples")'
        ),
        pms$devpars,
        outfile
    )
}


do_stats_perccellspersample = function(df_ncells, odir, pms) {
    outfile = file.path(odir, "perccellspersample.png")
    plotGG(
        df_ncells,
        "col",
        list(
            mapping = aes(
                x=seurat_clusters,
                y=CellFraction,
                fill=Sample
            ),
            position="stack"
        ),
        c(
            'ggtitle("Fraction of cells from each samples for each cluster")',
            'theme_prism()',
            'labs(x="Seurat Cluster", y="Fraction of total cells")'
        ),
        pms$devpars,
        outfile
    )
}


do_stats = function() {
    df_ncells = srtobj@meta.data |>
        group_by(Sample, seurat_clusters) |>
        summarize(nCells = n()) |>
        group_by(seurat_clusters) |>
        mutate(nCellsPerCluster = sum(nCells)) |>
        ungroup() |>
        group_by(Sample) |>
        mutate(
            percCellsPerSample = nCells / sum(nCells),
            nCellsPerSample = sum(nCells),
            CellFraction = nCells / nCellsPerCluster
        )

    # Sample seurat_clusters nCells nCellsPerCluster percCellsPerSample
    # <chr>  <fct>           <int>  <int>            <dbl>
    # s0301  0               397    2122             0.2296124928
    # s0301  1               211    1885             0.1220358589
    # ...
    #
    # nCellsPerSample CellFraction
    # <int>           <dbl>
    # 1729            397 / 2122
    # 1729            211 / 1885

    odir = file.path(outdir, "stats")
    dir.create(odir, showWarnings = FALSE)

    write.table(
        df_ncells,
        file = file.path(odir, "ncells.txt"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE
    )

    do_stats_ncells(df_ncells, odir, envs$stats$nCells)
    do_stats_ncellspersample(df_ncells, odir, envs$stats$nCellsPerSample)
    do_stats_perccellspersample(df_ncells, odir, envs$stats$percCellsPerSample)

}


do_exprs_ridgeplots = function(odir, pms) {
    outfile = file.path(odir, "ridgeplots.png")

    devpars = pms$devpars
    pms$devpars = NULL
    if (is.null(pms$ncol)) {
        pms$ncol = 2
    }
    if (is.null(pms$features)) {
        pms$features = VariableFeatures(srtobj)
    }
    if (is.null(devpars)) {
        devpars = list(
            width = 1000,
            height = ceiling(length(pms$features) / pms$ncol) * 500,
            res = 100
        )
    }
    pms$object = srtobj
    p = do.call(RidgePlot, pms)
    devpars$filename = outfile
    do.call(png, devpars)
    print(p)
    dev.off()
}


do_exprs_vlnplots = function(odir, pms) {
    outfile = file.path(odir, "vlnplots.png")

    devpars = pms$devpars
    pms$devpars = NULL
    boxplot = pms$boxplot
    pms$boxplot = NULL
    if (!is.null(boxplot) && length(boxplot) == 0) {
        boxplot = list(width = 0.1, fill = "white")
    }
    if (is.null(pms$ncol)) {
        pms$ncol = 2
    }
    if (is.null(pms$features)) {
        pms$features = VariableFeatures(srtobj)
    }
    if (is.null(devpars)) {
        devpars = list(
            width = 1000,
            height = ceiling(length(pms$features) / pms$ncol) * 500,
            res = 100
        )
    }
    pms$object = srtobj
    p = do.call(VlnPlot, pms)
    if (!is.null(boxplot)) {
        p = p + do.call(geom_boxplot, boxplot)
    }
    devpars$filename = outfile
    do.call(png, devpars)
    print(p)
    dev.off()
}


do_exprs_featureplots = function(odir, pms) {
    outfile = file.path(odir, "featureplots.png")

    devpars = pms$devpars
    pms$devpars = NULL
    if (is.null(pms$ncol)) {
        pms$ncol = 2
    }
    if (is.null(pms$features)) {
        pms$features = VariableFeatures(srtobj)
    }
    if (is.null(devpars)) {
        devpars = list(
            width = 1000,
            height = ceiling(length(pms$features) / pms$ncol) * 500,
            res = 100
        )
    }
    pms$object = srtobj
    p = do.call(FeaturePlot, pms)
    devpars$filename = outfile
    do.call(png, devpars)
    print(p)
    dev.off()
}

do_exprs_dotplot = function(odir, pms) {
    outfile = file.path(odir, "dotplot.png")

    devpars = pms$devpars
    pms$devpars = NULL
    plus = pms$plus
    pms$plus = NULL
    if (is.null(plus)) {
        plus = c()
    }
    if (is.null(pms$features)) {
        pms$features = VariableFeatures(srtobj)
    }
    if (is.null(devpars)) {
        devpars = list(
            height = length(unique(srtobj@meta.data$seurat_clusters)) * 80 + 150,
            width = length(pms$features) * 50 + 150,
            res = 100
        )
    }
    pms$object = srtobj
    p = do.call(DotPlot, pms)
    for (pls in plus) {
        p = p + eval(parse(text = pls))
    }
    devpars$filename = outfile
    do.call(png, devpars)
    print(p)
    dev.off()
}


do_exprs_heatmap = function(odir, pms) {
    outfile = file.path(odir, "heatmap.png")

    devpars = pms$devpars
    pms$devpars = NULL
    plus = pms$plus
    pms$plus = NULL
    if (is.null(plus)) {
        plus = c()
    }
    if (is.null(pms$features)) {
        pms$features = VariableFeatures(srtobj)
    }
    if (is.null(devpars)) {
        devpars = list(
            width = length(unique(srtobj@meta.data$seurat_clusters)) * 80 + 150,
            height = length(pms$features) * 50 + 150,
            res = 100
        )
    }
    downsample = pms$downsample
    pms$downsample = NULL
    if (is.null(downsample)) {
        pms$object = srtobj
    } else if (downsample %in% c("average", "mean")) {
        pms$object = AverageExpression(srtobj, return.seurat = TRUE)
    } else {
        pms$object = subset(srtobj, downsample = downsample)
    }
    p = do.call(DoHeatmap, pms)
    for (pls in plus) {
        p = p + eval(parse(text = pls))
    }
    devpars$filename = outfile
    do.call(png, devpars)
    print(p)
    dev.off()
}


do_exprs = function() {
    if (length(envs$exprs) == 0) {
        return (NULL)
    }
    odir = file.path(outdir, "exprs")
    dir.create(odir, showWarnings = FALSE)

    exprplots = names(envs$exprs)
    if ("ridgeplots" %in% exprplots) {
        do_exprs_ridgeplots(odir, envs$exprs$ridgeplots)
    }

    if ("vlnplots" %in% exprplots) {
        do_exprs_vlnplots(odir, envs$exprs$vlnplots)
    }

    if ("featureplots" %in% exprplots) {
        do_exprs_featureplots(odir, envs$exprs$featureplots)
    }

    if ("dotplot" %in% exprplots) {
        do_exprs_dotplot(odir, envs$exprs$dotplot)
    }

    if ("heatmap" %in% exprplots) {
        do_exprs_heatmap(odir, envs$exprs$heatmap)
    }
}


do_stats()
do_exprs()
