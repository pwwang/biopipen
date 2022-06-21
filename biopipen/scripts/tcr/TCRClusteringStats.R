source("{{biopipen_dir}}/utils/plot.R")
library(tibble)
library(tidyr)
library(dplyr)
library(immunarch)
library(ggprism)

immfile = {{in.immfile | quote}}
outdir = {{out.outdir | quote}}
envs = {{envs | r}}


cluster_size_distribution = function() {
    odir = file.path(outdir, "ClusterSizeDistribution")
    dir.create(odir, showWarnings = FALSE)

    clsizes = NULL
    for (sample in names(immdata$data)) {
        clsizes = bind_rows(
            clsizes,
            immdata$data[[sample]] |> count(TCR_Cluster) |> arrange(desc(n)) |> mutate(Sample = sample)
        )
    }

    outfile = file.path(odir, "cluster_size_distribution.txt")
    outplot = file.path(odir, "cluster_size_distribution.png")
    write.table(clsizes, outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

    plotGG(
        clsizes,
        "histogram",
        args = list(mapping = aes(x=n, fill=Sample)),
        ggs = c(
            "theme_prism()",
            "scale_y_continuous(trans='log10')",
            "labs(x='TCR cluster size', y='Count')"
        ),
        devpars = list(
            res = 100,
            height = 1000,
            width = 1000 + ceiling(length(immdata$data) / 16) * 150
        ),
        outfile = outplot
    )
}


shared_clusters = function() {
    odir = file.path(outdir, "SharedClusters")
    dir.create(odir, showWarnings = FALSE)
    tcr_clusters = list()
    samples = names(immdata$data)
    for (sample in samples) {
        tcr_clusters[[sample]] = immdata$data[[sample]] %>% pull("TCR_Cluster") %>% unique()
    }
    plotdata = matrix(NA, ncol = length(samples), nrow = length(samples))
    rownames(plotdata) = samples
    colnames(plotdata) = samples
    for (sample1 in samples) {
        for (sample2 in samples) {
            if (sample1 == sample2) {
                plotdata[sample1, sample2] = NA
            } else {
                plotdata[sample1, sample2] = length(
                    intersect(tcr_clusters[[sample1]], tcr_clusters[[sample2]])
                )
            }
        }
    }
    write.table(
        plotdata, file.path(odir, "shared_clusters.txt"),
        row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t"
    )

    if (is.null(envs$shared_clusters$heatmap_meta) || length(envs$shared_clusters$heatmap_meta) == 0) {
        anno = NULL
    } else {
        anno = as.list(immdata$meta[, heatmap_meta, drop=FALSE])
        anno = do.call(HeatmapAnnotation, anno)
    }

    # Plot heatmap
    plotHeatmap(
        plotdata,
        args = list(
            name = "Shared TCR Clusters",
            col = c("#ffe1e1", "red3"),
            cluster_rows = FALSE,
            top_annotation = anno,
            cell_fun = if (
                is.null(envs$shared_clusters$numbers_on_heatmap) || !envs$shared_clusters$numbers_on_heatmap
            ) NULL else function(j, i, x, y, width, height, fill) {
                grid.text(plotdata[samples[i], samples[j]], x, y, gp = gpar(fontsize = 10))
            }
        ),
        devpars = list(res=100, width=1000, height=1000),
        outfile = file.path(odir, "shared_clusters.png")

    )
}


shared_clusters_by_grouping = function() {
    if (is.null(envs$shared_clusters$grouping)) {
        return (NULL)
    }
    odir = file.path(outdir, "SharedClustersByGrouping")
    dir.create(odir, showWarnings = FALSE)

    data = list()
    grouping = envs$shared_clusters$grouping
    groups = immdata$meta |> pull(grouping) |> unique()
    sample_groups = list()
    for (group in groups) {
        for (sample in immdata$meta[
            immdata$meta[, grouping] == group,
            "Sample",
            drop = TRUE
        ]) {
            sample_groups[[sample]] = group
        }
        data[[group]] = c()
    }

    samples = names(immdata$data)
    for (sample in samples) {
        group = sample_groups[[sample]]
        clusters = immdata$data[[sample]] |> pull("TCR_Cluster") %>% unique()
        data[[group]] = unique(c(data[[group]], clusters))
    }

    outfile = file.path(odir, "shared_clusters_by_grouping.png")
    plotVenn(
        data,
        ggs = 'ggtitle("Shared TCR Clusters")',
        outfile = outfile
    )
}


sample_diversity = function() {
    odir = file.path(outdir, "SampleDiversity")
    dir.create(odir, showWarnings = FALSE)
    data = list()
    for (sample in names(immdata$data)) {
        data[[sample]] = immdata$data[[sample]] |> mutate(CDR3.aa = TCR_Cluster)
    }
    for (method in names(envs$sample_diversity)) {
        outfile = file.path(odir, paste0("diversity_", method, ".txt"))
        outplot = file.path(odir, paste0("diversity_", method, ".png"))
        div = repDiversity(data, .method = method)
        write.table(div, outfile, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
        if (method == "gini") {
            div = as.data.frame(div) |> rownames_to_column("Sample")
            colnames(div)[2] = "gini"
            div = left_join(div, immdata$meta, by="Sample")
            geom = "col"
            mapping = aes_string(
                x = "Sample",
                y = "gini",
                fill = "Sample"
            )
            ggs = c(
                "theme_prism(axis_text_angle = 90)",
                "labs(title='Gini coefficient', subtitle='Sample diversity estimation using the Gini coefficient')"
            )
            if (is.null(envs$sample_diversity[[method]]) || length(envs$sample_diversity[[method]]) == 0) {

            } else if (length(envs$sample_diversity[[method]]$by) == 1) {
                geom = "boxplot"
                mapping = aes_string(
                    x = envs$sample_diversity[[method]]$by,
                    y = "gini",
                    fill = envs$sample_diversity[[method]]$by
                )
            } else {
                div = div |> unite("Group", all_of(envs$sample_diversity[[method]]$by), sep="; ")
                geom = "boxplot"
                mapping = aes_string(
                    x = "Group",
                    y = "gini",
                    fill = "Group"
                )
            }
            plotGG(
                div,
                geom,
                args = list(mapping = mapping),
                ggs = ggs,
                outfile = outplot
            )

        } else {
            if (is.null(envs$sample_diversity[[method]]) || length(envs$sample_diversity[[method]]) == 0) {
                p = vis(div)
            } else {
                p = vis(
                    div,
                    .by = envs$sample_diversity[[method]]$by,
                    .meta = immdata$meta
                )
            }
            png(outplot, width=1000, height=1000, res=100)
            print(p)
            dev.off()
        }
    }
}


{
    # main
    # --------------------------------------------------
    # Load immunarch data
    immdata = readRDS(immfile)

    # Cluster size distribution
    cluster_size_distribution()

    # Shared clusters
    shared_clusters()

    # Shared clusters by grouping
    shared_clusters_by_grouping()

    # Diversity
    sample_diversity()
}
