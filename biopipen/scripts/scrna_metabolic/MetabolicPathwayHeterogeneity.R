source("{{biopipen_dir}}/utils/gsea.R")
source("{{biopipen_dir}}/utils/plot.R")

library(gtools)
library(parallel)
library(ggprism)

sceobjfile <- {{ in.sceobj | r }}
gmtfile <- {{ in.gmtfile | r }}
config <- {{ in.configfile | config: "toml" | r }}
outdir <- {{ out.outdir | r }}
envs <- {{envs | r}}

set.seed(8525)
groupby = config$grouping$groupby
if (grepl("^ident", groupby, ignore.case = TRUE)) {
    groupby = "seurat_clusters"
}

sceobj <- readRDS(sceobjfile)
do_one_subset <- function(subset) {
    subset_dir = file.path(outdir, subset)
    dir.create(subset_dir, showWarnings = FALSE)
    subset_sce <- sceobj[, sceobj$.subset == subset]
    metabolic_sce <- subset_sce[rowData(subset_sce)$metabolic, ]
    all_groups = as.character(metabolic_sce[[groupby]])
    groups <- unique(all_groups)

    enrich_data_df <- data.frame(x = NULL, y = NULL, NES = NULL, PVAL = NULL)
    pc_plotdata <- data.frame(
        x = numeric(), y = numeric(),
        sel = character(), group = character()
    )

    for (group in groups) {
        each_metabolic_sce <- metabolic_sce[, all_groups == group]
        each_metabolic_tpm <- assay(each_metabolic_sce, "exprs")
        each_metabolic_tpm <- each_metabolic_tpm[rowSums(each_metabolic_tpm) > 0, , drop=F]
        if (ncol(each_metabolic_tpm) == 1) {
            next
        }
        x <- each_metabolic_tpm
        ntop <- nrow(x)
        rv <- rowVars(x)
        select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        pca <- prcomp(t(x[select, ]))
        percentVar <- pca$sdev^2 / sum(pca$sdev^2)


        ### select PCs that explain at least 80% of the variance
        cum_var <- cumsum(percentVar)
        select_pcs <- which(cum_var > envs$select_pcs)[1]

        ### plot the PCA and explained variances
        tmp_plotdata <- data.frame(
            x = 1:length(percentVar), y = percentVar,
            sel = c(rep("y", select_pcs), rep("n", length(percentVar) - select_pcs)),
            group = rep(group, length(percentVar))
        )
        pc_plotdata <- rbind(pc_plotdata, tmp_plotdata)

        ####
        pre_rank_matrix <- as.matrix(rowSums(abs(pca$rotation[, 1:select_pcs, drop=FALSE])))
        pre_rank_matrix <- as.list(as.data.frame(t(pre_rank_matrix)))

        odir = file.path(subset_dir, group)
        dir.create(odir, showWarnings = FALSE)
        runFGSEA(
            pre_rank_matrix,
            gmtfile = gmtfile,
            top = 100,
            outdir = odir,
            plot = FALSE,
            envs = list(scoreType = "pos")
        )
        ############ Motify this
        result_file = file.path(odir, "fgsea.txt")
        gsea_result = read.table(result_file, header=T, row.names = NULL, sep="\t", check.names=F)
        # get the result
        enrich_data_df <- rbind(
            enrich_data_df,
            data.frame(x = group, y = gsea_result$pathway, NES = gsea_result$NES, PVAL = gsea_result$pval)
        )
    }

    # remove pvalue <0.01 pathways
    min_pval <- by(enrich_data_df$PVAL, enrich_data_df$y, FUN = min)
    select_pathways <- names(min_pval)[(min_pval <= envs$pathway_pval_cutoff)]
    select_enrich_data_df <- enrich_data_df[enrich_data_df$y %in% select_pathways, ]
    # converto pvalue to -log10
    pvals <- select_enrich_data_df$PVAL
    pvals[pvals <= 0] <- 1e-10
    select_enrich_data_df$PVAL <- -log10(pvals)

    # sort
    pathway_pv_sum <- by(select_enrich_data_df$PVAL, select_enrich_data_df$y, FUN = sum)
    pathway_order <- names(pathway_pv_sum)[order(pathway_pv_sum, decreasing = T)]
    ########################### top 10
    pathway_order <- pathway_order[1:10]
    select_enrich_data_df <- select_enrich_data_df[select_enrich_data_df$y %in% pathway_order, ]
    ########################################
    select_enrich_data_df$x <- factor(select_enrich_data_df$x, levels = mixedsort(groups))
    select_enrich_data_df$y <- factor(select_enrich_data_df$y, levels = pathway_order)

    ## buble plot
    select_enrich_data_df$x = sapply(select_enrich_data_df$x, function(x) {
        if (is.na(as.integer(x))) x else paste0("Cluster", x)
    })
    bubblefile = file.path(subset_dir, "pathway_heterogeneity.png")
    plotGG(
        select_enrich_data_df,
        "point",
        args = list(aes(x=x, y=y, size=PVAL, color=NES), shape=19),
        ggs = c(
            'scale_size(range = c(0, 5))',
            'scale_color_gradient(low = "white", high = "red")',
            'labs(
                x = NULL, y = NULL, color="NES", size="-log10(pval)"
            )',
            'theme_prism(axis_text_angle = 90)',
            'theme(legend.title = element_text())'
        ),
        devpars = envs$bubble_devpars,
        outfile = bubblefile
    )

    ## plot variance
    pc_plotdata$group <- factor(pc_plotdata$group, levels = mixedsort(groups))
    p <- ggplot(pc_plotdata) +
        geom_point(aes(x, y, colour = factor(sel)), size = 0.5) +
        scale_color_manual(values = c("gray", "#ff4000")) +
        facet_wrap(~group, scales = "free", ncol = 4) +
        theme_bw() +
        labs(x = "Principal components", y = "Explained variance (%)") +
        theme(
            legend.position = "none", panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(size = 0.2, colour = "black"),
            axis.ticks = element_line(colour = "black", size = 0.2),
            axis.text.x = element_text(colour = "black", size = 6),
            axis.text.y = element_text(colour = "black", size = 6),
            strip.background = element_rect(fill = "white", size = 0.2, colour = NULL),
            strip.text = element_text(size = 6)
        )

    ggsave(file.path(subset_dir, "PC_variance_plot.pdf"), p, device = "pdf", useDingbats = FALSE)
}

subsets <- unique(sceobj$.subset)

x = mclapply(subsets, do_one_subset, mc.cores = envs$ncores)
if (any(unlist(lapply(x, class)) == "try-error")) {
    stop("mclapply error")
}
