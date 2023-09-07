source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/plot.R")
library(tibble)
library(tidyr)
library(dplyr)
library(immunarch)
library(ggprism)

immfile = {{in.immfile | quote}}
outdir = {{out.outdir | quote}}
cluster_size_envs = {{envs.cluster_size | r}}
shared_clusters_envs = {{envs.shared_clusters | r}}
sample_diversity_envs = {{envs.sample_diversity | r}}

expand_cases = function(envs) {
    cases = envs$cases
    envs$cases = NULL
    if (is.null(cases) || length(cases) == 0) {
        cases = list(DEFAULT = list())
    }
    for (name in names(cases)) {
        case = cases[[name]]
        for (argname in names(envs)) {
            if (is.null(case[[argname]])) {
                case[[argname]] = envs[[argname]]
            }
        }
        for (n in c("height", "width", "res")) {
            if (is.null(case$devpars[[n]])) {
                case$devpars[[n]] = envs$devpars[[n]]
            }
        }
        if (is.null(case$devpars$height)) {
            case$devpars$height = 1000
        }
        if (is.null(case$devpars$width)) {
            case$devpars$width = 1000
        }
        if (is.null(case$devpars$res)) {
            case$devpars$res = 100
        }
        cases[[name]] = case
    }

    return (cases)
}

cluster_size_cases = expand_cases(cluster_size_envs)
shared_clusters_cases = expand_cases(shared_clusters_envs)
sample_diversity_cases = expand_cases(sample_diversity_envs)

cluster_size_distribution = function(name) {
    print(paste0("- Working on cluster size distribution: ", name))
    odir = file.path(outdir, "ClusterSizeDistribution", name)
    dir.create(odir, showWarnings = FALSE, recursive = TRUE)
    case = cluster_size_cases[[name]]

    clsizes = NULL
    for (sample in names(immdata$data)) {
        clsizes = bind_rows(
            clsizes,
            immdata$data[[sample]] %>% count(TCR_Cluster) %>% arrange(desc(n)) %>% mutate(Sample = sample)
        )
    }

    outfile = file.path(odir, "cluster_size_distribution.txt")
    outplot = file.path(odir, "cluster_size_distribution.png")
    write.table(clsizes, outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

    plotGG(
        clsizes,
        "histogram",
        args = list(mapping = aes_string(x="n", fill=case$by)),
        ggs = c(
            "theme_prism()",
            "scale_y_continuous(trans='log10')",
            "labs(x='TCR cluster size', y='Count')"
        ),
        devpars = case$devpars,
        outfile = outplot
    )
}

shared_clusters = function(name) {
    print(paste0("- Working on shared clusters: ", name))
    odir = file.path(outdir, "SharedClusters", name)
    dir.create(odir, showWarnings = FALSE, recursive = TRUE)
    case = shared_clusters_cases[[name]]
    if (!is.null(case$grouping)) {
        return(shared_clusters_by_grouping(name))
    }

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

    if (is.null(case$heatmap_meta) || length(case$heatmap_meta) == 0) {
        anno = NULL
    } else {
        anno = as.list(immdata$meta[, case$heatmap_meta, drop=FALSE])
        anno = do_call(ComplexHeatmap::HeatmapAnnotation, anno)
    }

    # Plot heatmap
    plotHeatmap(
        plotdata,
        args = list(
            name = "Shared TCR Clusters",
            col = c("#ffe1e1", "red3"),
            cluster_columns = FALSE,
            top_annotation = anno,
            cell_fun = if (
                is.null(case$numbers_on_heatmap) || !case$numbers_on_heatmap
            ) NULL else function(j, i, x, y, width, height, fill) {
                grid.text(plotdata[samples[i], samples[j]], x, y, gp = gpar(fontsize = 10))
            }
        ),
        devpars = case$devpars,
        outfile = file.path(odir, "shared_clusters.png")

    )
}

shared_clusters_by_grouping = function(name) {
    odir = file.path(outdir, "SharedClusters", name)
    case = shared_clusters_cases[[name]]

    data = list()
    grouping = case$grouping
    groups = immdata$meta %>% pull(grouping) %>% unique()
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
        clusters = immdata$data[[sample]] %>% pull("TCR_Cluster") %>% unique()
        data[[group]] = unique(c(data[[group]], clusters))
    }

    outfile = file.path(odir, "shared_clusters.png")
    plotVenn(
        data,
        ggs = 'ggtitle("Shared TCR Clusters")',
        devpars = case$devpars,
        outfile = outfile
    )
}


sample_diversity = function(name) {
    print(paste0("- Working on sample diversity: ", name))
    odir = file.path(outdir, "SampleDiversity", name)
    dir.create(odir, showWarnings = FALSE, recursive = TRUE)
    case = sample_diversity_cases[[name]]

    data = list()
    for (sample in names(immdata$data)) {
        data[[sample]] = immdata$data[[sample]] %>% mutate(CDR3.aa = TCR_Cluster)
    }
    outfile = file.path(odir, "diversity.txt")
    outplot = file.path(odir, "diversity.png")
    div = repDiversity(data, .method = case$method)
    write.table(div, outfile, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
    if (case$method == "gini") {
        div = as.data.frame(div) %>% rownames_to_column("Sample")
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
        if (is.null(case$by) || length(case$by) == 0) {

        } else {
            case$by = trimws(strsplit(case$by, ",")[[1]])
            if (length(case$by) == 1) {
                geom = "boxplot"
                mapping = aes_string(x = case$by, y = "gini", fill = case$by)
            } else {
                div = div %>% unite("Group", all_of(case$by), sep="; ")
                geom = "boxplot"
                mapping = aes_string(
                    x = "Group",
                    y = "gini",
                    fill = "Group"
                )
            }
        }

        plotGG(
            div,
            geom,
            args = list(mapping = mapping),
            ggs = ggs,
            devpars = case$devpars,
            outfile = outplot
        )

    } else {
        if (is.null(case$by) || length(case$by) == 0) {
            p = vis(div)
        } else {
            p = vis(
                div,
                .by = trimws(strsplit(case$by, ",")[[1]]),
                .meta = immdata$meta
            )
        }
        png(
            outplot,
            width=case$devpars$width, height=case$devpars$height, res=case$devpars$res
        )
        print(p)
        dev.off()
    }
}


{
    # main
    # --------------------------------------------------
    # Load immunarch data
    immdata = readRDS(immfile)

    # Cluster size distribution
    sapply(names(cluster_size_cases), cluster_size_distribution)

    # Shared clusters
    sapply(names(shared_clusters_cases), shared_clusters)

    # Diversity
    sapply(names(sample_diversity_cases), sample_diversity)
}
