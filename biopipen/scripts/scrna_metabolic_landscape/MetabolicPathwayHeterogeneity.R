source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/gsea.R")
source("{{biopipen_dir}}/utils/plot.R")

library(gtools)
library(parallel)
library(ggprism)
library(Matrix)
library(sparseMatrixStats)
library(Seurat)

sobjfile <- {{ in.sobjfile | r }}
outdir <- {{ out.outdir | r }}
joboutdir <- {{ job.outdir | r }}
gmtfile <- {{ envs.gmtfile | r }}
select_pcs <- {{ envs.select_pcs | r }}
ncores <- {{ envs.ncores | r }}
pathway_pval_cutoff <- {{ envs.pathway_pval_cutoff | r }}
bubble_devpars <- {{ envs.bubble_devpars | r }}
grouping <- {{ envs.grouping | r }}
grouping_prefix <- {{ envs.grouping_prefix | r }}
subsetting_cols <- {{ envs.subsetting | r }}
subsetting_prefix <- {{ envs.subsetting_prefix | r }}

if (!is.null(grouping_prefix) && nchar(grouping_prefix) > 0) {
    grouping_prefix = paste0(grouping_prefix, "_")
}

if (!is.null(subsetting_prefix) && nchar(subsetting_prefix) > 0) {
    subsetting_prefix = paste0(subsetting_prefix, "_")
}

set.seed(8525)

## gmt_pathways is copied from fgsea package.
gmt_pathways <- function(gmt_file) {
    pathway_lines <- strsplit(readLines(gmt_file), "\t")
    pathways <- lapply(pathway_lines, tail, -2)
    names(pathways) <- sapply(pathway_lines, head, 1)
    pathways
}

gmtfile <- localizeGmtfile(gmtfile)
pathways <- gmt_pathways(gmtfile)
metabolics <- unique(as.vector(unname(unlist(pathways))))
sobj <- readRDS(sobjfile)


do_one_subset <- function(s, subset_col, subset_prefix) {
    log_info(paste0("  Handling subset value: ", s, " ..."))
    if (is.null(s)) {
        subset_dir = file.path(outdir, "ALL")
        subset_obj = sobj
    } else {
        subset_dir = file.path(outdir, slugify(paste0(subset_prefix, s)))
        subset_code = paste0("subset(sobj, subset = ", subset_col, " == '", s, "')")
        subset_obj = eval(parse(text = subset_code))
    }
    dir.create(subset_dir, showWarnings = FALSE)

    features = intersect(rownames(subset_obj), metabolics)
    all_groups = as.character(subset_obj@meta.data[[grouping]])
    groups <- unique(all_groups)

    enrich_data_df <- data.frame(x = NULL, y = NULL, NES = NULL, PVAL = NULL)
    pc_plotdata <- data.frame(
        x = numeric(),
        y = numeric(),
        sel = character(),
        group = character()
    )

    for (group in groups) {
        group_code = paste0("subset(subset_obj, subset = ", grouping, " == '", group, "')")
        each_metabolic_obj <- eval(parse(text = group_code))
        each_metabolic_exprs <- GetAssayData(each_metabolic_obj)[features, , drop=F]
        each_metabolic_exprs <- each_metabolic_exprs[rowSums(each_metabolic_exprs) > 0, , drop=F]
        if (ncol(each_metabolic_exprs) == 1) { next }
        x <- each_metabolic_exprs
        ntop <- nrow(x)
        rv <- rowVars(x)
        select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        pca <- prcomp(Matrix::t(x[select, ]))
        percentVar <- pca$sdev^2 / sum(pca$sdev^2)

        ### select PCs that explain at least 80% of the variance
        cum_var <- cumsum(percentVar)
        selected_pcs <- which(cum_var > select_pcs)[1]

        ### plot the PCA and explained variances
        tmp_plotdata <- data.frame(
            x = 1:length(percentVar),
            y = percentVar,
            sel = c(rep("y", selected_pcs), rep("n", length(percentVar) - selected_pcs)),
            group = rep(group, length(percentVar))
        )
        pc_plotdata <- rbind(pc_plotdata, tmp_plotdata)

        ####
        pre_rank_matrix <- as.matrix(rowSums(abs(pca$rotation[, 1:selected_pcs, drop=FALSE])))
        pre_rank_matrix <- as.list(as.data.frame(t(pre_rank_matrix)))

        odir = file.path(subset_dir, paste0(grouping_prefix, slugify(group)))
        dir.create(odir, showWarnings = FALSE)
        runFGSEA(
            pre_rank_matrix,
            gmtfile = gmtfile,
            top = 100,
            outdir = odir,
            plot = FALSE,
            envs = list(scoreType = "std", nproc=1)
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

    # remove pvalue < 0.01 pathways
    min_pval <- by(enrich_data_df$PVAL, enrich_data_df$y, FUN = min)
    select_pathways <- names(min_pval)[(min_pval <= pathway_pval_cutoff)]
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
    select_enrich_data_df$x = sapply(select_enrich_data_df$x, function(x) { paste0(grouping_prefix, x) })
    bubblefile = file.path(subset_dir, "pathway_heterogeneity.png")
    bub_devpars = list() # bubble_devpars
    if (is.null(bub_devpars$res)) {
        bub_devpars$res = 100
    }
    if (is.null(bub_devpars$width)) {
        bub_devpars$width = 300 +
            max(nchar(as.character(select_enrich_data_df$y))) * 8 +
            length(unique(select_enrich_data_df$x)) * 25
    }
    if (is.null(bub_devpars$height)) {
        bub_devpars$height = 400 +
            max(nchar(unique(select_enrich_data_df$x))) * 8 +
            length(unique(select_enrich_data_df$y)) * 25
    }
    bub_devpars$height = max(bub_devpars$height, 480)
    # For debug purposes
    write.table(
        select_enrich_data_df,
        file.path(subset_dir, "pathway_heterogeneity.txt"),
        sep="\t",
        quote=F,
        row.names=F
    )
    if (nrow(select_enrich_data_df) == 0) {
        p = ggplot(data.frame(text = "No significant pathways found")) +
            geom_text(aes(x = 0, y = 0, label = text), size = 10) +
            theme_void() +
            theme(
                plot.margin = unit(c(0, 0, 0, 0), "cm"),
                plot.background = element_rect(fill = "white", colour = NA)
            )
        png(bubblefile, width = 600, height = 100, res = 70)
        print(p)
        dev.off()
    } else {
        plotGG(
            select_enrich_data_df,
            "point",
            args = list(aes(x=x, y=y, size=PVAL, color=NES), shape=19),
            ggs = c(
                'scale_size(range = c(2, 10))',
                'scale_color_gradient(low = "white", high = "red")',
                'labs(
                    x = NULL, y = NULL, color="NES", size="-log10(pval)"
                )',
                'theme_prism(axis_text_angle = 90)',
                'theme(legend.title = element_text())'
            ),
            devpars = bub_devpars,
            outfile = bubblefile
        )
    }

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
            axis.line = element_line(linewidth = 0.2, colour = "black"),
            axis.ticks = element_line(colour = "black", linewidth = 0.2),
            axis.text.x = element_text(colour = "black", size = 6),
            axis.text.y = element_text(colour = "black", size = 6),
            strip.background = element_rect(fill = "white", linewidth = 0.2, colour = NULL),
            strip.text = element_text(size = 6)
        )

    ggsave(file.path(subset_dir, "PC_variance_plot.pdf"), p, device = "pdf", useDingbats = FALSE)

    list(
        list(kind = "descr", content = "Metabolic pathways enriched in genes with highest contribution to the metabolic heterogeneities"),
        list(kind = "image", src = bubblefile),
        h1 = ifelse(is.null(s), "Metabolic pathway heterogeneity", paste0(subset_prefix, s))
    )
}

do_one_subset_col <- function(subset_col, subset_prefix) {
    log_info(paste0("- Handling subset column: ", subset_col, " ..."))
    if (is.null(subset_col)) {
        x <- do_one_subset(NULL, subset_col = NULL, subset_prefix = NULL)
        do.call(add_report, x)
    } else {
        subsets <- na.omit(unique(sobj@meta.data[[subset_col]]))

        if (ncores == 1) {
            x = lapply(subsets, do_one_subset, subset_col = subset_col, subset_prefix = subset_prefix)
        } else {
            x <- mclapply(subsets, do_one_subset, subset_col = subset_col, subset_prefix = subset_prefix, mc.cores = ncores)
            if (any(unlist(lapply(x, class)) == "try-error")) {
                stop(paste0("\nmclapply error:", x))
            }
        }
        for (r in x) { do.call(add_report, r) }
    }
}

if (is.null(subsetting_cols)) {
    do_one_subset_col(NULL)
} else {
    for (i in seq_along(subsetting_cols)) {
        do_one_subset_col(subsetting_cols[i], subsetting_prefix[i])
    }
}

save_report(joboutdir)
