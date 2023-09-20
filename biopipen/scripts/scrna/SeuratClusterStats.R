source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/plot.R")
library(Seurat)
library(rlang)
library(dplyr)
library(tibble)
library(ggprism)
library(tidyseurat)

srtfile = {{in.srtobj | r}}
outdir = {{out.outdir | r}}
envs = {{envs | r: todot="-"}}

srtobj = readRDS(srtfile)

do_stats_cells = function(casename, devpars, odir, by = NULL, frac = FALSE, filtering = NULL) {
    plotfile = file.path(odir, paste0(casename, ".png"))
    txtfile = paste0(plotfile, ".txt")
    if (is.null(devpars)) {
        devpars = list(res = 100, height = 1000)
    } else {
        devpars$res = if (is.null(devpars$res)) default_devpars$res else devpars$res
        devpars$height = if (is.null(devpars$height)) default_devpars$height else devpars$height
    }
    if (frac) {
        ylab = paste("Fraction of cells")
        mapping_y = "cellFraction"
    } else {
        ylab = paste("Number of cells")
        mapping_y = "nCells"
    }
    df_cells = srtobj@meta.data
    if (!is.null(filtering)) {
        df_cells = df_cells %>% filter(!!rlang::parse_expr(filtering))
    }
    if (is.null(by)) {
        df_cells = df_cells %>%
            group_by(seurat_clusters) %>%
            summarize(nCells = n(), .groups = "keep") %>%
            mutate(cellFraction = nCells / sum(nCells))

        if (is.null(devpars$width)) { devpars$width = 1000 }
        plotGG(
            df_cells,
            geom = "col",
            args = list(
                mapping = aes_string(x="seurat_clusters", y=mapping_y)
            ),
            ggs = c(
                paste0('ggtitle("', ylab, ' for each cluster")'),
                'theme_prism()',
                'theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))',
                paste0('labs(x="Seurat Cluster", y="', ylab, '")')
            ),
            devpars = devpars,
            outfile = plotfile
        )
    } else {
        df_cells = df_cells %>%
            group_by(!!sym(by), seurat_clusters) %>%
            summarize(nCells = n()) %>%
            group_by(seurat_clusters) %>%
            mutate(cellFraction = nCells / sum(nCells))

        if (is.null(devpars$width)) {
            devpars$width = 800 + 200 * ceiling(length(unique(df_cells[[by]])) / 20)
        }
        plotGG(
            df_cells,
            geom = "col",
            args = list(
                mapping = aes_string(x="seurat_clusters", y=mapping_y, fill=by),
                position = "stack"
            ),
            ggs = c(
                paste0('ggtitle("', ylab, ' for each cluster by ', by, '")'),
                'theme_prism()',
                'theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))',
                paste0('labs(x="Seurat Cluster", y="', ylab, '")')
            ),
            devpars = devpars,
            outfile = plotfile
        )
    }

    write.table(
        df_cells,
        file = txtfile,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE
    )
}

do_stats = function() {
    odir = file.path(outdir, "stats")
    dir.create(odir, showWarnings = FALSE)
    for (name in names(envs$stats)) {
        stat_pars = envs$stats[[name]]
        args = list(
            casename = name,
            devpars = stat_pars$devpars,
            odir = odir,
            by = stat_pars$by,
            frac = FALSE,
            filtering = stat_pars$filter
        )
        if (startsWith(name, "fracCells")) {
            args$frac = TRUE
        } else if (!startsWith(name, "nCells")) {
            warning(paste("Unknown stat:", name, ", skipping"))
            next
        }

        do_call(do_stats_cells, args)
    }
}

.get_outfile = function(odir, prefix, ext = "png") {
    i = 1
    while (TRUE) {
        outfile = file.path(odir, paste0(prefix, "-", i, ".", ext))
        if (!file.exists(outfile)) {
            return(outfile)
        }
        i = i + 1
    }
    return(outfile)
}

.get_features = function(features, genes, default = NULL) {
    if (is.null(default)) {
        default = VariableFeatures(srtobj)[1:20]
    }
    # When nothing passed, use the genes
    if (is.null(features)) {
        if (is.null(genes)) {
            return (default)
        } else {
            return (genes)
        }
    }
    # When multiple items passed, use them as features
    if (length(features) > 1) {
        return (features)
    }
    # See if it is "default"
    if (features == "default") {
        return (default)
    }
    # See if it is a file
    if (!file.exists(features)) {
        return (features)
    }
    # length(features) == 1 && file.exists(features[1])
    feats = read.table(features, header = FALSE, row.names = NULL, check.names = FALSE)
    feats$V1
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
    pms$features = .get_features(pms$features, genes)
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
        pms$object = srtobj %>% filter(eval(parse(text=subsetpms)))
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
        boxplot = list(width = .1, fill = "white")
    }
    if (is.null(pms$ncol)) {
        pms$ncol = 2
    }
    pms$features = .get_features(pms$features, genes)
    if (is.null(devpars)) {
        devpars = list(
            width = pms$ncol * 480,
            height = ceiling(length(pms$features) / pms$ncol) * 480,
            res = 100
        )
    }
    if (is.null(subsetpms)) {
        pms$object = srtobj
    } else {
        pms$object = srtobj %>% tidyseurat::filter(eval(parse(text=subsetpms)))
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
    pms$features = .get_features(pms$features, genes)
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
        pms$object = srtobj %>% filter(eval(parse(text=subsetpms)))
    }
    p = do_call(FeaturePlot, pms)
    devpars$filename = outfile

    tryCatch({
        do_call(png, devpars)
        print(p)
        dev.off()
    }, error = function(e) {
        stop(
            paste(
                paste(names(devpars), collapse=" "),
                paste(devpars, collapse=" "),
                e,
                sep = "\n"
            )
        )
    })
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
    pms$features = .get_features(pms$features, genes)
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
        pms$object = srtobj %>% filter(eval(parse(text=subsetpms)))
    }
    p = do_call(DotPlot, pms)
    for (pls in plus) {
        p = p + eval(parse(text = pls))
    }
    devpars$filename = outfile
    tryCatch({
        do_call(png, devpars)
        print(p)
        dev.off()
    }, error = function(e) {
        stop(
            paste(
                paste(names(devpars), collapse=" "),
                paste(devpars, collapse=" "),
                e,
                sep = "\n"
            )
        )
    })
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
    pms$features = .get_features(pms$features, genes)
    if (is.null(devpars)) {
        devpars = list(
            width = length(unique(srtobj@meta.data$seurat_clusters)) * 60 + 150,
            height = length(pms$features) * 40 + 150,
            res = 100
        )
    }
    downsample = pms$downsample
    pms$downsample = NULL

    if (is.null(subsetpms)) {
        sobj = srtobj
    } else {
        sobj = srtobj %>% filter(eval(parse(text=subsetpms)))
    }
    if (is.null(downsample)) {
        pms$object = sobj
        warn(
            paste0(
                "DoHeatmap: `downsample` not specified, using full data. ",
                "This may cause a blank heatmap. ",
                "See: https://github.com/satijalab/seurat/issues/2724"
            ),
            immediate. = TRUE
        )
    } else if (downsample %in% c("average", "mean")) {
        pms$object = AverageExpression(sobj, return.seurat = TRUE)
    } else {
        pms$object = subset(sobj, downsample = downsample)
    }
    p = do_call(DoHeatmap, pms)
    for (pls in plus) {
        p = p + eval(parse(text = pls))
    }
    mapal = colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
    p = p + scale_fill_gradientn(colours = rev(mapal))
    devpars$filename = outfile
    do_call(png, devpars)
    print(p)
    dev.off()
}


do_exprs_table = function(odir, pms, genes) {
    outfile = .get_outfile(odir, "table", "tsv")

    subsetpms = pms$subset
    log2_scale = pms$log2
    if (is.null(log2_scale)) { log2_scale = TRUE }
    features = .get_features(pms$features, genes)
    title = pms$title
    if (is.null(title)) { title = tools::file_path_sans_ext(basename(outfile)) }
    cat(title, file = paste0(outfile, ".title"))

    if (is.null(subsetpms)) {
        sobj = srtobj
    } else {
        sobj = srtobj %>% filter(eval(parse(text=subsetpms)))
    }
    # default slot (data), assay
    # values are exponentiated prior to averaging so that averaging is done in non-log space.
    avgdata = AverageExpression(sobj, features = features)
    edata = avgdata$RNA
    # replace the missing genes
    edata[rownames(avgdata$integrated), ] = avgdata$integrated
    if (log2_scale) { edata = log2(edata) }
    edata = as.data.frame(edata) %>%
        rownames_to_column("Gene") %>%
        select(Gene, everything())
    write.table(edata, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
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
        genes = read.table(genes, header = FALSE, sep = "\t", row.names = NULL, check.names = FALSE)
        genes = genes[,1,drop=TRUE]
    } else if (!is.null(genes) && is.character(genes)) {
        genes = trimws(strsplit(genes, ",")[[1]])
    }

    exprplots = names(envs$exprs)
    for (name in exprplots) {
        cat(paste0("Expr plot: ", name, " ...\n"), file = stderr())
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
        } else if (startsWith(name, "table")) {
            do_exprs_table(odir, envs$exprs[[name]], genes)
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
