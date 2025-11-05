library(rlang)
library(glue)
library(dplyr)
library(scplotter)
library(biopipen.utils)

screpfile <- {{in.screpfile | r}}
outdir <- {{out.outdir | r}}
joboutdir <- {{job.outdir | r}}
envs <- {{envs | r: todot="-"}}
mutaters <- envs$mutaters
cases <- envs$cases
envs$mutaters <- NULL
envs$cases <- NULL

log <- get_logger()
reporter <- get_reporter()

VIZ_TYPE_TO_SECTION <- list(
    volume = "Number of Clones",
    abundance = "Clonal Abundance",
    length = "Clonal Sequence Length",
    residency = "Clonal Residency",
    stat = "Clonal Statistics",
    composition = "Clonal Composition",
    overlap = "Clonal Overlap",
    diversity = "Clonal Diversity",
    geneusage = "Gene Usage",
    positional = "Positional Properties",
    kmer = "Kmer Analysis",
    rarefaction = "Rarefaction Analysis"
)

get_plot_descr <- function(viz_type, case) {
    if (identical(viz_type, "volume")) {
        if (identical(case$plot_type %||% "bar", "bar")) {
            out <- glue(
                "This bar graph illustrates the distribution of unique clones across {x}(s). ",
                "The x-axis represents the different {x}(s), while the y-axis denotes the number of unique clones.",
                x = case$x %||% "Sample")
        } else {  # box/violin
            out <- glue(
                "This {case$plot_type} plot compares the distribution of unique clones across {x}(s). ",
                "The x-axis represents {x}(s), and the data points were broken down by samples, while the y-axis denotes the number of unique clones.",
                x = case$x)
            if (!is.null(case$comparisons)) {
                out <- glue("{out} The p-values of the comparisons are performed using {case$pairwise_method %||% 'wilcox.test'}.")
            }
        }
        out <- glue("{out} The clones are identified by {case$clone_call %||% 'aa'} and {case$chain %||% 'both'} chain(s) was/were used.")
    } else if (identical(viz_type, "abundance")) {
        if ((case$plot_type %||% "trend") %in% c("trend", "histogram")) {
            out <- glue(
                "The abundance plot illustrates the number of clones at different abundance levels. ",
                "The x-axis represents the abundance levels (different sizes of clones), while the y-axis denotes the number of clones at each level."
            )
        } else {  # density
            out <- glue(
                "The abundance plot illustrates the distribution of clones at different abundance levels. ",
                "The x-axis represents the abundance levels (different sizes of clones), while the y-axis denotes the density of clones at each level."
            )
        }
        out <- glue("{out} The clones are identified by {case$clone_call %||% 'aa'} and {case$chain %||% 'both'} chain(s) was/were used.")
    } else if (identical(viz_type, "length")) {
        if (identical(case$plot_type %||% "bar", "bar")) {
            out <- glue(
                "This bar graph illustrates the distribution of the length of the CDR3 sequences. ",
                "The x-axis represents the different lengths of the CDR3 sequences, while the y-axis denotes the number of corresponding sequences. ",
                "{case$chain %||% 'both'} chain(s) was/were counted."
            )
        } else if (identical(case$plot_type, "density")) {
            out <- glue(
                "This density plot illustrates the distribution of the length of the CDR3 sequences. ",
                "The x-axis represents the different lengths of the CDR3 sequences, while the y-axis denotes the density of corresponding sequences. ",
                "{case$chain %||% 'both'} chain(s) was/were counted."
            )
        } else {  # box/violin
            out <- glue(
                "This {case$plot_type} plot compares the distribution of the length of the CDR3 sequences. ",
                "The x-axis represents the different lengths of the CDR3 sequences, while the y-axis denotes the number of corresponding sequences. ",
                "{case$chain %||% 'both'} chain(s) was/were counted."
            )
            if (!is.null(case$comparisons)) {
                out <- glue("{out} The p-values of the comparisons are performed using {case$pairwise_method %||% 'wilcox.test'}.")
            }
        }
    } else if (identical(viz_type, "residency")) {
        if (identical(case$plot_type %||% "scatter", "scatter")) {
            out <- glue(
                "This scatter plot illustrates the residency of clones across different groups (axes). ",
                "The x-axis represents the first group, while the y-axis denotes the second group. ",
                "The size of the points represents the size of the clones. For the shared clones, the size of the points ",
                "represents the {case$scatter_size_by %||% 'max'} size of the clones in the groups. ",
                "The color of the points represents the category of the clones. 'Singlet' are the clones with a single cell; ",
                "'Expanded' are the clones with multiple cells; 'Dual' are the clones with cells in both groups; ",
                "'Dual (g1 > g2)' are the clones with cells in both groups and more cells in the first group; ",
                "'Dual (g1 < g2)' are the clones with cells in both groups and more cells in the second group; and ",
                "'Dual (Equal)' are the clones with cells in both groups and equal cells in both groups. ",
                "The {case$scatter_cor %||% 'pearson'} correlation was calculated based on the size of the shared clones and p-value ",
                "is also shown in the subtitle."
            )
        } else if (identical(case$plot_type, "venn")) {
            out <- glue(
                "This {case$plot_type} plot illustrates the residency of clones across different groups. ",
                "The categories of the clones are shown in the Venn diagram, and the color represents the size the category. ",
                "The number of singlets are also annotated in the plot."
            )
        } else {  # upset
            out <- glue(
                "This {case$plot_type} plot illustrates the residency of clones across different groups. ",
                "The categories of the clones are shown as rows in the bottom table of the plot. The intersections of the categories ",
                "are shown as the connected lines in the table and the size of the intersections are shown as the bars on the top of ",
                "the plot. The color of the bars represents the size of the intersections."
            )
        }
        out <- glue("{out} The clones are identified by {case$clone_call %||% 'aa'} and {case$chain %||% 'both'} chain(s) was/were used.")
    } else if (identical(viz_type, "stat")) {
        if (case$plot_type %in% c("sankey", "alluvial")) {
            out <- glue(
                "This {case$plot_type} plot illustrates the statistics of clones across different groups. ",
                "The bars are showing the groups and the flow/links are showing the transitions of the clones. "
            )
        } else {  # trend
            out <- glue(
                "This trend plot illustrates the statistics of clones across different groups. ",
                "The x-axis represents the groups, while the y-axis denotes the number/fraction of clones. ",
                "The links between the groups are showing the transitions of the clones. "
            )
        }
        out <- glue("{out} The clones are identified by {case$clone_call %||% 'aa'} and {case$chain %||% 'both'} chain(s) was/were used.")
    } else if (identical(viz_type, "composition")) {
        plot_type <- case$plot_type %||% "bar"
        method <- case$method %||% "homeostasis"
        if (plot_type %in% c("bar", "ring")) {
            out <- glue("This {plot_type} graph illustrates the composition of the categories of the clones.")
        } else {  # box/violin
            out <- glue(
                "This {plot_type} plot compares the composition of the categories of the clones. ",
                "For each category in the x-axis, the values are from each sample. "
            )
            if (!is.null(case$comparisons)) {
                out <- glue("{out} The p-values of the comparisons are performed using {case$pairwise_method %||% 'wilcox.test'}.")
            }
        }
        if (method %in% c("homeostasis", "homeo", "rel")) {
            out <- glue("{out} The clone categories are defined based on the relative abundance of the clones in the samples.")
        } else if (method == "top") {
            out <- glue(
                "{out} The clone categories are defined based on the size of the clones in the samples. ",
                "The largest clone ranks as 1, and clones are marked by their indexes."
            )
        } else if (method == "rare") {
            out <- glue(
                "{out} The clone categories are defined based on the size of the clones in the samples. ",
                "The clones are categorized literally based on the size of the clones. For example, ",
                "1 means the singlet clones."
            )
        }
        out <- glue("{out} The clones are identified by {case$clone_call %||% 'aa'} and {case$chain %||% 'both'} chain(s) was/were used.")
    } else if (identical(viz_type, "overlap")) {
        out <- glue(
            "This heatmaps illustrates the overlap of clones across different groups, using \"{case$method %||% 'raw'}\" as the metric. "
        )
        if ((case$method %||% "raw") == "raw") {
            out <- glue("{out} The 'raw' metric is the raw size of the intersection of the clones. ")
        } else if (case$method == "overlap") {
            out <- glue(
                "{out} The 'overlap' metric, a.k.a. 'overlap coefficient' or 'Szymkiewicz-Simpson coefficient', ",
                "is calculated as the size of the intersection divided by the size of the smaller set. "
            )
        } else if (case$method == "morisita") {
            out <- glue("{out} The 'morisita' metric, a.k.a. 'Morisita`s overlap index', ",
                "named after Masaaki Morisita, is a statistical measure of dispersion of individuals in a population. ",
                "This formula is based on the assumption that increasing the size of the samples will increase the diversity ",
                "because it will include different habitats. The value is 0 if the two samples do not overlap in terms of species, ",
                "and 1 if the species occur in the same proportions in both samples."
            )
        } else if (case$method == "jaccard") {
            out <- glue(
                "{out} The 'jaccard' metric, a.k.a. 'Jaccard index', ",
                "is a statistic used for gauging the similarity and diversity of sample sets. ",
                "It is defined in general taking the ratio of two sizes (areas or volumes), ",
                "the intersection size divided by the union size, also called intersection over union."
            )
        } else if (case$method == "cosine") {
            out <- glue("{out} The 'cosine' metric, a.k.a. 'cosine similarity', ",
                "is calculated as the dot product of the vectors divided by the product of the magnitudes of the vectors. "
            )
        }
        out <- glue("{out} The clones are identified by {case$clone_call %||% 'aa'} and {case$chain %||% 'both'} chain(s) was/were used.")
    } else if (identical(viz_type, "diversity")) {
        if (identical(case$plot_type %||% "bar", "bar")) {
            out <- glue(
                "This bar graph illustrates the diversity of the clones across different groups. ",
                "The x-axis represents the groups, while the y-axis denotes the diversity of the clones. "
            )
        } else {  # box/violin
            out <- glue(
                "This {case$plot_type} plot compares the diversity of the clones across different groups. ",
                "The x-axis represents the groups, while the y-axis denotes the diversity of the clones. ",
                "For each group on the x-axis, the values are calculated from each sample. "
            )
            if (!is.null(case$comparisons)) {
                out <- glue("{out} The p-values of the comparisons are performed using {case$pairwise_method %||% 'wilcox.test'}.")
            }
        }
        if (identical(case$method %||% "shannon", "shannon")) {
            out <- glue(
                "{out} The diversity is calculated using the Shannon index, ",
                "a statistic that estimates the diversity of a biological community by considering both species richness and evenness. ",
                "The index is calculated as the negative sum of the product of the proportion of each species and the ",
                "logarithm of that proportion. If practically all abundance is concentrated to one type, and the other ",
                "types are very rare (even if there are many of them), Shannon entropy approaches zero. ",
                "When there is only one type in the dataset, Shannon entropy exactly equals zero."
            )
        } else if (identical(case$method, "inv.simpson")) {
            out <- glue(
                "{out} The diversity is calculated using the Inverse Simpson index, simply equals true diversity of order 2. ",
                "It is calculated as the reciprocal of the sum of the squares of the proportion of each species. "
            )
        } else if (identical(case$method, "norm.entropy")) {
            out <- glue(
                "{out} The diversity is calculated using the Normalized Shannon entropy, ",
                "a statistic that estimates the diversity of a biological community by considering both species richness and evenness. ",
                "The index is calculated as the negative sum of the product of the proportion of each species and the ",
                "logarithm of that proportion. If practically all abundance is concentrated to one type, and the other ",
                "types are very rare (even if there are many of them), Shannon entropy approaches zero. ",
                "When there is only one type in the dataset, Shannon entropy exactly equals zero. ",
                "The normalized Shannon entropy is calculated as the Shannon entropy divided by the maximum possible entropy. ",
                "The maximum possible entropy is the Shannon entropy when all species are equally abundant."
            )
        } else if (identical(case$method, "gini.simpson")) {
            out <- glue(
                "{out} The diversity is calculated using the Gini-Simpson index, ",
                "a.k.a Gini impurity. The original Simpson index λ equals the probability that two entities taken at random ",
                "from the dataset of interest (with replacement) represent the same type. Its transformation 1 - λ, ",
                "therefore, equals the probability that the two entities represent different types."
            )
        } else if (identical(case$method, "chao1")) {
            out <- glue(
                "{out} The diversity is calculated using the Chao1 index. ",
                "The Chao1 index is a lower bound estimate of the true species richness of a community. ",
                "It is based on the number of singletons and doubletons in a sample. ",
                "The Chao1 index is calculated as the observed number of species plus the square of the number of singlets ",
                "divided by twice the number of doublets. "
            )
        } else if (identical(case$method, "ACE")) {
            out <- glue(
                "{out} The diversity is calculated using the Abundance-based Coverage Estimator (ACE) index. ",
                "a statistical method used to estimate the number of species in a community or ecosystem. ",
                "The ACE index is calculated as the sum of the estimated number of species in the sample. ",
                "The ACE index is based on the abundance of species and estimates the number of additional ",
                "species that are likely to be present but have not been observed. ",
                "Higher ACE values indicate higher species diversity."
            )
        } else if (identical(case$method, "gini.coeff")) {
            out <- glue(
                "{out} The diversity is calculated using the Gini coefficient. ",
                "The Gini coefficient is a measure of statistical dispersion intended to represent the income or wealth distribution of a nation's residents. ",
                "The Gini coefficient ranges from 0 to 1, where 0 represents perfect equality and 1 represents perfect inequality. ",
                "A Gini coefficient of 0 means that all individuals have the same income, while a Gini coefficient of 1 means that one individual has all the income and the rest have none."
            )
        } else if (case$method %in% c("d50", "dXX")) {
            out <- glue(
                "{out} The diversity is calculated using the {case$method} index. ",
                "The {case$method} index a metric used to measure the diversity of an immune repertoire by calculating the ",
                "percentage of unique clonotype that account for the top {case[['d']] %||% 50}% of all clones. "
            )
        }
        out <- glue("{out} The clones are identified by {case$clone_call %||% 'aa'} and {case$chain %||% 'both'} chain(s) was/were used.")
    } else if (viz_type == "geneusage") {
        if (identical(case$plot_type %||% "bar", "bar")) {
            out <- glue(
                "This bar graph illustrates the distribution of the gene usage of the clones. ",
                "The x-axis represents the different genes, while the y-axis denotes the fraction of the genes. "
            )
        } else if (identical(case$plot_type, "heatmap")) {
            out <- glue(
                "This heatmap illustrates the distribution of the gene usage of the clones in different groups. ",
                "The color of the heatmap represents the fraction of the genes. "
            )
        } else if (case$plot_type %in% c("circos", "chord")) {
            out <- glue("This {case$plot_type} plot illustrates the distribution of the gene usage of the clones")
            if (length(case$genes) == 1) {
                out <- glue(
                    "{out} for the '{case$genes[[1]]}' genes in different groups. ",
                    "The links are showing the composition of the genes in the groups. "
                )
            } else {
                out <- glue(
                    "{out} for the {paste(case$genes, collapse = ' and ')} genes. ",
                    "The links are showing the connections between the genes. "
                )
            }
        } else {  # sankey/alluvial
            out <- glue(
                "This {case$plot_type} plot illustrates the distribution of the gene usage of the clones in different groups. ",
                "The bars are showing the groups and the flow/links are showing the transitions of the gene usages. "
            )
        }
    } else if (viz_type == "positional") {
        if ((case$plot_type %||% "bar") %in% c("bar", "line")) {
            out <- glue(
                "This {case$plot_type %||% 'bar'} graph illustrates the distribution of the positional properties of the amino acids in the CDR3 sequences. "
            )
        } else if (identical(case$plot_type, "heatmap")) {
            out <- glue(
                "This heatmap illustrates the distribution of the positional properties of the amino acids in the CDR3 sequences ",
                "in different groups. The color of the heatmap represents the intensity of the properties. "
            )
        } else if (case$plot_type %in% c("violin", "box")) {
            out <- glue(
                "This {case$plot_type} plot compares the distribution of the positional properties of the amino acids in the CDR3 sequences. ",
                "The x-axis represents the different positions of the amino acids, while the y-axis denotes the intensity of the properties. ",
                "For each position on the x-axis, the values are calculated from each sample. "
            )
            if (!is.null(case$comparisons)) {
                out <- glue("{out} The p-values of the comparisons are performed using {case$pairwise_method %||% 'wilcox.test'}.")
            }
        }
        if (identical(case$method %||% "AA", "AA")) {
            out <- glue("{out} The positional properties are calculated based on amino acid frequency at each position.")
        } else if (identical(case$method, "shannon")) {
            out <- glue(
                "{out} The positional properties are calculated based on the Shannon entropy of the amino acids at each position. ",
                "The Shannon entropy is a statistic that estimates the diversity of a biological community by considering both species richness and evenness. ",
                "The index is calculated as the negative sum of the product of the proportion of each species and the ",
                "logarithm of that proportion. If practically all abundance is concentrated to one type, and the other ",
                "types are very rare (even if there are many of them), Shannon entropy approaches zero."
            )
        } else if (identical(case$method, "inv.simpson")) {
            out <- glue(
                "{out} The positional properties are calculated based on the Inverse Simpson index of the amino acids at each position. ",
                "The Inverse Simpson index is calculated as the reciprocal of the sum of the squares of the proportion of each species. "
            )
        } else if (identical(case$method, "norm.entropy")) {
            out <- glue(
                "{out} The positional properties are calculated based on the Normalized Shannon entropy of the amino acids at each position. ",
                "The Normalized Shannon entropy is a statistic that estimates the diversity of a biological community by considering both species richness and evenness. ",
                "The index is calculated as the negative sum of the product of the proportion of each species and the ",
                "logarithm of that proportion. If practically all abundance is concentrated to one type, and the other ",
                "types are very rare (even if there are many of them), Shannon entropy approaches zero. ",
                "When there is only one type in the dataset, Shannon entropy exactly equals zero. ",
                "The normalized Shannon entropy is calculated as the Shannon entropy divided by the maximum possible entropy. ",
                "The maximum possible entropy is the Shannon entropy when all species are equally abundant."
            )
        } else if (identical(case$method, "Atchley")) {
            out <- glue(
                "{out} The positional properties are calculated based on the Atchley factors of the amino acids at each position. ",
                "The Atchley factors are a set of six factors that describe the amino acids based on their physicochemical properties. ",
                "The factors are polarity, secondary structure, volume, codon diversity, electrostatic charge, and solvent accessibility. "
            )
        } else if (identical(case$method, "Kidera")) {
            out <- glue(
                "{out} The positional properties are calculated based on the Kidera factors of the amino acids at each position. ",
                "The Kidera factors are a set of ten factors that describe the amino acids based on their physicochemical properties. ",
                "The factors are hydrophobicity, hydrophilicity, side chain mass, pK1, pK2, pI, side chain pK, alpha helix, beta sheet, and turn. "
            )
        } else if (identical(case$method, "stScales")) {
            out <- glue(
                "{out} The positional properties are calculated based on the stScales of the amino acids at each position. ",
                "The stScales were proposed by Yang et al, taking 827 properties into account which are mainly constitutional, ",
                "topological, geometrical, hydrophobic, elec- tronic, and steric properties. "
            )
        } else if (identical(case$method, "tScales")) {
            out <- glue(
                "{out} The positional properties are calculated based on the tScales of the amino acids at each position. ",
                "The tScales are a descriptor set for amino acids that are based on topological descriptors. "
            )
        } else if (identical(case$method, "VHSE")) {
            out <- glue(
                "{out} The positional properties are calculated based on Vectors of Hydrophobic, Steric, and Electronic (VHSE) properties of the amino acids at each position. "
            )
        }
        out <- glue("{out} The properties are calculated based on the {case$chain %||% 'both'} chain(s).")
    } else if (viz_type == "kmer") {
        if ((case$plot_type %||% "bar") %in% c("bar", "line")) {
            out <- glue(
                "This {case$plot_type %||% 'bar'} graph illustrates the distribution of the kmers of the CDR3 sequences. "
            )
        } else if (case$plot_type == "heatmap") {
            out <- glue(
                "This heatmap illustrates the distribution of the kmers of the CDR3 sequences in different groups. ",
                "The color of the heatmap represents the frequency of the kmers. "
            )
        }
        out <- glue("{out} The kmers are calculated based on the {case$k %||% 3}-mers of the CDR3 sequences, and {case$chain %||% 'both'} chain(s) was/were used.")
    } else if (viz_type == "rarefaction") {
        out <- glue(
            "This rarefaction curve illustrates the number of unique clones as a function of the number of cells. ",
            "The x-axis represents the number of cells, while the y-axis denotes the number of unique clones. ",
            "The solid line represents the observed number of unique clones, while the dashed line represents the expected number of unique clones. ",
            "The clones are identified by {case$clone_call %||% 'aa'} and {case$chain %||% 'both'} chain(s) was/were used."
        )
    }

    out
}

log$info("Loading scRepertoire object ...")
screp <- read_obj(screpfile)

log$info("Applying mutaters if any ...")
screp <- ScRepMutate(screp, mutaters)

log$info("Making cases ...")
cases <- expand_cases(cases, envs)
viz_types <- list()
for (name in names(cases)) {
    case <- cases[[name]]
    if (is.null(case$viz_type)) {
        stop("Error: Visualization type is not defined for case '", name, "'")
    }
    if (!case$viz_type %in% names(VIZ_TYPE_TO_SECTION)) {
        stop("Error: Unknown visualization type '", case$viz_type, "' for case '", name,
             "'. Available types: ", paste(names(VIZ_TYPE_TO_SECTION), collapse = ", "))
    }
    if (!grepl("::", name, fixed = TRUE)) {
        viz_types[[case$viz_type]] <- viz_types[[case$viz_type]] %||% 0
        viz_types[[case$viz_type]] <- viz_types[[case$viz_type]] + 1
    }
}
cases <- list_rename(cases, function(name, case) {
    if (!grepl("::", name, fixed = TRUE) && viz_types[[case$viz_type]] > 1) {
        section <- VIZ_TYPE_TO_SECTION[[case$viz_type]]
        return(paste0(section, "::", name))
    }
    return(TRUE)
})

do_case <- function(name, case) {
    log$info("- Case: {name}")
    info <- case_info(name, outdir, is_dir = FALSE, create = TRUE)

    case <- extract_vars(case, "viz_type", "descr", "devpars", "more_formats", "save_code", subset_ = "subset")

    if (!is.null(subset_)) {
        case$data <- ScRepSubset(screp, subset_)
    } else {
        case$data <- screp
    }
    fnname <- tools::toTitleCase(viz_type)
    if (fnname == "Geneusage") {
        fnname <- "GeneUsage"
    }
    plot_fn <- paste0("Clonal", fnname, "Plot")
    plot_fn <- utils::getFromNamespace(plot_fn, "scplotter")
    if (is.null(plot_fn)) {
        stop("Error: Unknown visualization type: ", viz_type)
    }

    p <- do_call(plot_fn, case)
    save_plot(p, info$prefix, devpars, formats = unique(c("png", more_formats)))

    report <- list(
        kind = "table_image",
        src = paste0(info$prefix, ".png"),
        download = list(),
        descr = html_escape(descr %||% get_plot_descr(viz_type, case)),
        name = html_escape(info$name)
    )
    exformats <- setdiff(more_formats, "png")
    if (length(exformats) > 0) {
        report$download <- lapply(exformats, function(fmt) {
            paste0(info$prefix, ".", fmt)
        })
    }

    if (isTRUE(save_code)) {
        save_plotcode(
            p,
            setup = c('library(scplotter)', '', 'load("data.RData")'),
            prefix = info$prefix,
            "case"
        )
        report$download <- c(report$download, list(list(
            src = paste0(info$prefix, ".code.zip"),
            tip = "Download the code to reproduce the plot",
            icon = "Code"
        )))
    }

    reporter$add2(report, hs = c(info$section, info$name), ui = "table_of_images:2")
}

lapply(names(cases), function(name) do_case(name, cases[[name]]))

reporter$save(joboutdir)
