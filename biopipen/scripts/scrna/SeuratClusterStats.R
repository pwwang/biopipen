source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/plot.R")
library(Seurat)
library(dplyr)
library(ggprism)
library(tidyseurat)

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

.get_outfile = function(odir, prefix) {
    i = 1
    while (TRUE) {
        outfile = file.path(odir, paste0(prefix, "-", i, ".png"))
        if (!file.exists(outfile)) {
            return(outfile)
        }
        i = i + 1
    }
    return(outfile)
}


do_exprs_ridgeplots = function(odir, pms, genes) {
    outfile = .get_outfile(odir, "ridgeplots")

    devpars = pms$devpars
    pms$devpars = NULL
    if (is.null(pms$ncol)) {
        pms$ncol = 2
    }
    plus = pms$plus
    pms$plus = NULL
    title = pms$title
    pms$title = NULL
    if (is.null(title)) {
        title = tools::file_path_sans_ext(basename(outfile))
    }
    cat(title, file = paste0(outfile, ".title"))
    subsetpms = pms$subset
    pms$subset = NULL
    if (is.null(plus)) {
        plus = c()
    }
    if (is.null(pms$features) && is.null(genes)) {
        pms$features = VariableFeatures(srtobj)[1:20]
    } else if (is.null(pms$features) && !is.null(genes)) {
        pms$features = genes
    }
    if (is.null(devpars)) {
        devpars = list(
            width = 1000,
            height = ceiling(length(pms$features) / pms$ncol) * 250,
            res = 100
        )
    }
    if (is.null(subsetpms)) {
        pms$object = srtobj
    } else {
        pms$object = srtobj |> filter(eval(parse(text=subsetpms)))
    }
    p = do_call(RidgePlot, pms)
    for (pls in plus) {
        p = p + eval(parse(text = pls))
    }
    devpars$filename = outfile
    do_call(png, devpars)
    print(p)
    dev.off()
}


do_exprs_vlnplots = function(odir, pms, genes) {
    outfile = .get_outfile(odir, "vlnplots")

    devpars = pms$devpars
    pms$devpars = NULL
    boxplot = pms$boxplot
    pms$boxplot = NULL
    plus = pms$plus
    pms$plus = NULL
    subsetpms = pms$subset
    pms$subset = NULL
    title = pms$title
    pms$title = NULL
    if (is.null(title)) {
        title = tools::file_path_sans_ext(basename(outfile))
    }
    cat(title, file = paste0(outfile, ".title"))
    if (is.null(plus)) {
        plus = c()
    }
    if (!is.null(boxplot) && length(boxplot) == 0) {
        boxplot = list(width = 0.1, fill = "white")
    }
    if (is.null(pms$ncol)) {
        pms$ncol = 2
    }
    if (is.null(pms$features) && is.null(genes)) {
        pms$features = VariableFeatures(srtobj)[1:20]
    } else if (is.null(pms$features) && !is.null(genes)) {
        pms$features = genes
    }
    if (is.null(devpars)) {
        devpars = list(
            width = 1000,
            height = ceiling(length(pms$features) / pms$ncol) * 250,
            res = 100
        )
    }
    if (is.null(subsetpms)) {
        pms$object = srtobj
    } else {
        pms$object = srtobj |> filter(eval(parse(text=subsetpms)))
    }
    p = do_call(VlnPlot, pms)
    if (!is.null(boxplot)) {
        p = p + do_call(geom_boxplot, boxplot)
    }
    for (pls in plus) {
        p = p + eval(parse(text = pls))
    }
    devpars$filename = outfile
    do_call(png, devpars)
    print(p)
    dev.off()
}


do_exprs_featureplots = function(odir, pms, genes) {
    outfile = .get_outfile(odir, "featureplots")

    devpars = pms$devpars
    pms$devpars = NULL
    title = pms$title
    pms$title = NULL
    if (is.null(title)) {
        title = tools::file_path_sans_ext(basename(outfile))
    }
    cat(title, file = paste0(outfile, ".title"))
    if (is.null(pms$ncol)) {
        pms$ncol = 2
    }
    if (is.null(pms$features) && is.null(genes)) {
        pms$features = VariableFeatures(srtobj)[1:20]
    } else if (is.null(pms$features) && !is.null(genes)) {
        pms$features = genes
    }
    if (is.null(devpars)) {
        devpars = list(
            width = 1000,
            height = ceiling(length(pms$features) / pms$ncol) * 250,
            res = 100
        )
    }
    subsetpms = pms$subset
    pms$subset = NULL
    if (is.null(subsetpms)) {
        pms$object = srtobj
    } else {
        pms$object = srtobj |> filter(eval(parse(text=subsetpms)))
    }
    p = do_call(FeaturePlot, pms)
    devpars$filename = outfile
    do_call(png, devpars)
    print(p)
    dev.off()
}

do_exprs_dotplot = function(odir, pms, genes) {
    outfile = .get_outfile(odir, "dotplot")

    devpars = pms$devpars
    pms$devpars = NULL
    plus = pms$plus
    pms$plus = NULL
    subsetpms = pms$subset
    pms$subset = NULL
    title = pms$title
    pms$title = NULL
    if (is.null(title)) {
        title = tools::file_path_sans_ext(basename(outfile))
    }
    cat(title, file = paste0(outfile, ".title"))
    if (is.null(plus)) {
        plus = c()
    }
    if (is.null(pms$features) && is.null(genes)) {
        pms$features = VariableFeatures(srtobj)[1:20]
    } else if (is.null(pms$features) && !is.null(genes)) {
        pms$features = genes
    }
    if (is.null(devpars)) {
        devpars = list(
            height = length(unique(srtobj@meta.data$seurat_clusters)) * 80 + 150,
            width = length(pms$features) * 50 + 150,
            res = 100
        )
    }
    if (is.null(subsetpms)) {
        pms$object = srtobj
    } else {
        pms$object = srtobj |> filter(eval(parse(text=subsetpms)))
    }
    p = do_call(DotPlot, pms)
    for (pls in plus) {
        p = p + eval(parse(text = pls))
    }
    devpars$filename = outfile
    do_call(png, devpars)
    print(p)
    dev.off()
}


do_exprs_heatmap = function(odir, pms, genes) {
    outfile = .get_outfile(odir, "heatmap")

    devpars = pms$devpars
    pms$devpars = NULL
    plus = pms$plus
    pms$plus = NULL
    subsetpms = pms$subset
    pms$subset = NULL
    title = pms$title
    pms$title = NULL
    if (is.null(title)) {
        title = tools::file_path_sans_ext(basename(outfile))
    }
    cat(title, file = paste0(outfile, ".title"))
    if (is.null(plus)) {
        plus = c()
    }
    if (is.null(pms$features) && is.null(genes)) {
        pms$features = VariableFeatures(srtobj)[1:20]
    } else if (is.null(pms$features) && !is.null(genes)) {
        pms$features = genes
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

    if (is.null(subsetpms)) {
        sobj = srtobj
    } else {
        sobj = srtobj |> filter(eval(parse(text=subsetpms)))
    }
    if (is.null(downsample)) {
        pms$object = sobj
    } else if (downsample %in% c("average", "mean")) {
        pms$object = AverageExpression(sobj, return.seurat = TRUE)
    } else {
        pms$object = subset(sobj, downsample = downsample)
    }
    p = do_call(DoHeatmap, pms)
    for (pls in plus) {
        p = p + eval(parse(text = pls))
    }
    devpars$filename = outfile
    do_call(png, devpars)
    print(p)
    dev.off()
}


do_exprs = function() {
    if (length(envs$exprs) == 0) {
        return (NULL)
    }
    odir = file.path(outdir, "exprs")
    dir.create(odir, showWarnings = FALSE)
    genes = envs$exprs$genes
    envs$exprs$genes = NULL
    if (!is.null(genes) && is.character(genes) && file.exists(genes)) {
        genes = read.table(genes, header = FALSE, row.names = NULL, check.names = FALSE)
        genes = genes[,1,drop=TRUE]
    }

    exprplots = names(envs$exprs)
    for (name in exprplots) {
        if (startsWith(name, "ridgeplots")) {
            do_exprs_ridgeplots(odir, envs$exprs[[name]], genes)
        } else if (startsWith(name, "vlnplots")) {
            do_exprs_vlnplots(odir, envs$exprs[[name]], genes)
        } else if (startsWith(name, "featureplots")) {
            do_exprs_featureplots(odir, envs$exprs[[name]], genes)
        } else if (startsWith(name, "dotplot")) {
            do_exprs_dotplot(odir, envs$exprs[[name]], genes)
        } else if (startsWith(name, "heatmap")) {
            do_exprs_heatmap(odir, envs$exprs[[name]], genes)
        } else {
            print(paste("Unrecognized expression plot type: ", name))
        }
    }
}

do_dimplot = function(odir, dpname, dpenvs) {
    devpars = dpenvs$devpars
    dpenvs$devpars = NULL
    if (is.null(devpars)) {
        devpars = list(
            width = 1000,
            height = 1000,
            res = 100
        )
    }
    plus = dpenvs$plus
    dpenvs$plus = NULL
    if (is.null(plus)) {
        plus = c()
    }
    if (!any(grepl("ggtitle(", plus, fixed=T))) {
        plus = c(plus, paste0("ggtitle(", dQuote(dpname, q=F), ")"))
    }

    dpenvs$object = srtobj
    p = do_call(DimPlot, dpenvs)
    for (pls in plus) {
        p = p + eval(parse(text = pls))
    }

    devpars$filename = file.path(odir, paste0(slugify(dpname), ".png"))
    do_call(png, devpars)
    print(p)
    dev.off()
}

do_dimplots = function() {
    if (length(envs$dimplots) == 0) {
        return (NULL)
    }
    odir = file.path(outdir, "dimplots")
    dir.create(odir, showWarnings = FALSE)

    for (dpname in names(envs$dimplots)) {
        dpenvs = envs$dimplots[[dpname]]
        do_dimplot(odir, dpname, dpenvs)
    }
}


do_stats()
do_exprs()
do_dimplots()
