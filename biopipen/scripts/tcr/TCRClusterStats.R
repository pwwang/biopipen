{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "plot.R" | source_r }}
library(tibble)
library(tidyr)
library(dplyr)
library(rlang)
library(immunarch)
library(ggprism)

immfile = {{in.immfile | r}}
outdir = {{out.outdir | r}}
cluster_size_envs = {{envs.cluster_size | r}}
shared_clusters_envs = {{envs.shared_clusters | r}}
sample_diversity_envs = {{envs.sample_diversity | r}}
joboutdir = {{job.outdir | r}}

log_info("Expanding analysis cases ...")
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
    log_info("- Working on cluster size distribution: {name}")

    odir = file.path(outdir, "ClusterSizeDistribution", slugify(name))
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
    outplot_pdf = file.path(odir, "cluster_size_distribution.pdf")
    write.table(clsizes, outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

    plotGG(
        clsizes,
        "histogram",
        args = list(mapping = aes(x=n, fill=!!sym(case$by))),
        ggs = c(
            "theme_prism()",
            "scale_y_continuous(trans='log10')",
            "labs(x='TCR cluster size', y='Count')",
            "scale_fill_biopipen()"
        ),
        devpars = case$devpars,
        outfile = c(outplot, outplot_pdf)
    )

    add_report(
        list(
            src = outplot,
            name = ifelse(name == "DEFAULT", FALSE, name),
            descr = paste0("Cluster size distribution for each ", case$by),
            download = outplot_pdf
        ),
        ui = "table_of_images",
        h1 = "Cluster Size Distribution"
    )
}

shared_clusters = function(name) {
    log_info("- Working on shared clusters: {name}")

    odir = file.path(outdir, "SharedClusters", slugify(name))
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

    if (!is.null(case$sample_order) && length(case$sample_order) > 0) {
        if (length(case$sample_order) == 1) {
            case$sample_order = trimws(strsplit(case$sample_order, ",")[[1]])
        }
        nonexisting = setdiff(case$sample_order, samples)
        if (length(nonexisting) > 0) {
            stop(paste("  The following samples do not exist in `sample_order`:", paste(nonexisting, collapse=", ")))
        }
        plotdata = plotdata[, case$sample_order, drop=FALSE]
    }

    if (is.null(case$heatmap_meta) || length(case$heatmap_meta) == 0) {
        anno = NULL
    } else {
        anno = as.list(
            immdata$meta[
                match(colnames(plotdata), immdata$meta$Sample),
                case$heatmap_meta,
                drop=FALSE
            ])
        anno = do_call(ComplexHeatmap::HeatmapAnnotation, anno)
    }

    cluster_rows = case$cluster_rows && nrow(plotdata) > 2
    col_samples = colnames(plotdata)
    if (!cluster_rows) {
        plotdata = plotdata[col_samples, ]
        row_samples = col_samples
    } else {
        row_samples = samples
    }

    hmplot = file.path(odir, "shared_clusters.png")
    hmplot_pdf = file.path(odir, "shared_clusters.pdf")
    # Plot heatmap
    plotHeatmap(
        plotdata,
        args = list(
            name = "Shared TCR Clusters",
            col = c("#ffe1e1", "red3"),
            cluster_columns = FALSE,
            cluster_rows = cluster_rows,
            top_annotation = anno,
            cell_fun = if (
                is.null(case$numbers_on_heatmap) || !case$numbers_on_heatmap
            ) NULL else function(j, i, x, y, width, height, fill) {
                grid.text(row_samples[i], col_samples[j], x, y, gp = gpar(fontsize = 10))
            }
        ),
        devpars = case$devpars,
        outfile = c(hmplot, hmplot_pdf)
    )

    add_report(
        list(
            src = hmplot,
            download = hmplot_pdf,
            name = ifelse(name == "DEFAULT", FALSE, name),
            descr = paste0("Shared TCR clusters across samples")
        ),
        ui = "table_of_images",
        h1 = "Shared TCR Clusters"
    )
}

shared_clusters_by_grouping = function(name) {
    odir = file.path(outdir, "SharedClusters", slugify(name))
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
    outfile_pdf = file.path(odir, "shared_clusters.pdf")
    plotVenn(
        data,
        ggs = 'ggtitle("Shared TCR Clusters")',
        devpars = case$devpars,
        outfile = c(outfile, outfile_pdf)
    )

    add_report(
        list(
            src = outfile,
            download = outfile_pdf,
            name = ifelse(name == "DEFAULT", FALSE, name),
            descr = paste0("Shared TCR clusters across ", grouping)
        ),
        ui = "table_of_images",
        h1 = "Shared TCR Clusters"
    )
}


div_methods = list(
    gini = list(
        name = "The Gini coefficient",
        descr = "The Gini coefficient is a measure of statistical dispersion intended to represent the income or wealth distribution of a nation's residents, and is the most commonly used measurement of inequality."
    ),
    gini.simp = list(
        name = "The Gini-Simpson index",
        descr = "The Gini-Simpson index is a measure of diversity. It is one of the most commonly used in ecology. It is also known as the Simpson index, the Simpson concentration index, the Simpson dominance index, or the Simpson diversity index."
    ),
    inv.simp = list(
        name = "The inverse Simpson index",
        descr = "It is the effective number of types that is obtained when
                 the weighted arithmetic mean is used to quantify average
                 proportional abundance of types in the dataset of interest."
    ),
    div = list(
        name = "The true diversity",
        descr = "It refers to the number of equally abundant types needed
                 for the average proportional abundance of the types to
                 equal that observed in the dataset of interest where all
                 types may not be equally abundant."
    )
)

sample_diversity = function(name) {
    log_info("- Working on sample diversity: {name}")

    odir = file.path(outdir, "SampleDiversity", slugify(name))
    dir.create(odir, showWarnings = FALSE, recursive = TRUE)
    case = sample_diversity_cases[[name]]

    data = list()
    for (sample in names(immdata$data)) {
        data[[sample]] = immdata$data[[sample]] %>% mutate(CDR3.aa = TCR_Cluster)
    }
    outfile = file.path(odir, "diversity.txt")
    outplot = file.path(odir, "diversity.png")
    outplot_pdf = file.path(odir, "diversity.pdf")
    div = repDiversity(data, .method = case$method)
    write.table(
        if (ncol(div) == 1) {
            as.data.frame(div) %>% rownames_to_column("Sample")
        } else {
            div
        },
        outfile,
        row.names=TRUE,
        col.names=TRUE,
        quote=FALSE,
        sep="\t"
    )

    if (case$method == "gini") {
        div = as.data.frame(div) %>% rownames_to_column("Sample")
        colnames(div)[2] = "gini"
        div = left_join(div, immdata$meta, by="Sample")
        geom = "col"
        mapping = aes(x = Sample, y = gini, fill = Sample)
        ggs = c(
            "theme_prism(axis_text_angle = 90)",
            "labs(title='Gini coefficient', subtitle='Sample diversity estimation using the Gini coefficient')",
            "scale_fill_biopipen()"
        )
        if (is.null(case$by) || length(case$by) == 0) {

        } else {
            case$by = trimws(strsplit(case$by, ",")[[1]])
            if (length(case$by) == 1) {
                geom = "boxplot"
                mapping = aes(x = !!sym(case$by), y = gini, fill = !!sym(case$by))
            } else {
                div = div %>% unite("Group", all_of(case$by), sep="; ")
                geom = "boxplot"
                mapping = aes(x = Group, y = gini, fill = Group)
            }
        }

        plotGG(
            div,
            geom,
            args = list(mapping = mapping),
            ggs = ggs,
            devpars = case$devpars,
            outfile = c(outplot, outplot_pdf)
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

        pdf(
            outplot_pdf,
            width=case$devpars$width / case$devpars$res,
            height=case$devpars$height / case$devpars$res
        )
        print(p)
        dev.off()
    }

    add_report(
        list(
            ui = "flat",
            label = "Diversity Plot",
            contents = list(
                list(
                    kind = "descr",
                    content = paste(
                        div_methods[[case$method]]$name,
                        ifelse(
                            is.null(case$by) || length(case$by) == 0,
                            "",
                            paste0(" grouped by ", paste(case$by, collapse = ", "))
                        ),
                        div_methods[[case$method]]$descr
                    )
                ),
                list(
                    kind = "image",
                    src = outplot,
                    download = outplot_pdf
                )
            )
        ),
        list(
            ui = "flat",
            label = "Diversity Table",
            contents = list(
                list(kind = "table", src = outfile, data = list(index_col = 0))
            )
        ),
        ui = "tabs",
        h2 = ifelse(name == "DEFAULT", "#", name),
        h1 = "Sample Diversity using TCR clusters"
    )
}


{
    # main
    # --------------------------------------------------
    # Load immunarch data
    log_info("Loading immunarch data ...")
    immdata = readRDS(immfile)

    # Cluster size distribution
    log_info("Performing cluster size distribution analysis ...")
    sapply(names(cluster_size_cases), cluster_size_distribution)

    # Shared clusters
    log_info("Performing shared clusters analysis ...")
    sapply(names(shared_clusters_cases), shared_clusters)

    # Diversity
    log_info("Performing sample diversity analysis ...")
    sapply(names(sample_diversity_cases), sample_diversity)

    save_report(joboutdir)
}
