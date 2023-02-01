source("{{biopipen_dir}}/utils/plot.R")

library(scater)
library(reshape2)
library(RColorBrewer)
library(parallel)
library(ggprism)
library(Seurat)

sobjfile <- {{ in.sobjfile | r }}
outdir <- {{ out.outdir | r }}
gmtfile <- {{ envs.gmtfile | r }}
ntimes <- {{ envs.ntimes | r }}
ncores <- {{ envs.ncores | r }}
heatmap_devpars <- {{ envs.heatmap_devpars | r }}
violin_devpars <- {{ envs.violin_devpars | r }}
grouping <- {{ envs.grouping | r }}
grouping_prefix <- {{ envs.grouping_prefix | r }}
subsetting <- {{ envs.subsetting | r }}
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

pathways <- gmt_pathways(gmtfile)
pathway_names <- names(pathways)
metabolics <- unique(as.vector(unname(unlist(pathways))))
sobj <- readRDS(sobjfile)
DefaultAssay(sobj) <- "RNA"

## calculate how many pathways of one gene involved.
num_of_pathways <- function(gmtfile, overlapgenes) {
    filter_pathways <- list()
    for (p in pathway_names) {
        genes <- pathways[[p]]
        common_genes <- intersect(genes, overlapgenes)
        if (length(common_genes >= 5)) {
            filter_pathways[[p]] <- common_genes
        }
    }

    all_genes <- unique(as.vector(unlist(filter_pathways)))
    gene_times <- data.frame(num = rep(0, length(all_genes)), row.names = all_genes)
    for (p in pathway_names) {
        for (g in filter_pathways[[p]]) {
            gene_times[g, "num"] <- gene_times[g, "num"] + 1
        }
    }
    gene_times
}

do_one_subset <- function(s) {
    if (is.null(s)) {
        subset_dir <- file.path(outdir, "ALL")
        dir.create(subset_dir, showWarnings = FALSE)
        subset_obj <- sobj
    } else {
        subset_dir <- file.path(outdir, paste0(subsetting_prefix, s))
        dir.create(subset_dir, showWarnings = FALSE)

        subset_code = paste0("subset(sobj, subset = ", subsetting, " == '", s, "')")
        subset_obj = eval(parse(text = subset_code))
    }

    all_cell_types <- subset_obj@meta.data[[grouping]]
    cell_types <- unique(all_cell_types)

    gene_pathway_number <- num_of_pathways(
        gmtfile,
        intersect(rownames(subset_obj), metabolics)
    )

    ## Calculate the pathway activities
    # mean ratio of genes in each pathway for each cell type
    mean_expression_shuffle <- matrix(
        NA,
        nrow = length(pathway_names),
        ncol = length(cell_types),
        dimnames = list(pathway_names, cell_types)
    )
    mean_expression_noshuffle <- matrix(
        NA,
        nrow = length(pathway_names),
        ncol = length(cell_types),
        dimnames = list(pathway_names, cell_types)
    )
    ### calculate the pvalues using shuffle method
    pvalues_mat <- matrix(
        NA,
        nrow = length(pathway_names),
        ncol = length(cell_types),
        dimnames = (list(pathway_names, cell_types))
    )

    for (p in pathway_names) {
        genes <- pathways[[p]]
        genes_comm <- intersect(genes, rownames(subset_obj))
        genes_expressed <- names(rowSums(subset_obj)[rowSums(subset_obj) > 0])
        genes_comm <- intersect(genes_comm, genes_expressed)
        if (length(genes_comm) < 5) next

        pathway_metabolic_obj <- subset(subset_obj, features = genes_comm)
        mean_exp_eachCellType <- AverageExpression(pathway_metabolic_obj)$RNA

        # remove genes which are zeros in any celltype to avoid extreme ratio value
        keep <- rownames(mean_exp_eachCellType)[rowAlls(mean_exp_eachCellType > 0.001)]
        if (length(keep) < 3) next

        # using the loweset value to replace zeros for avoiding extreme ratio value
        pathway_metabolic_obj <- subset(pathway_metabolic_obj, features = keep)
        assay_data = GetAssayData(pathway_metabolic_obj)
        assay_data <- t(apply(assay_data, 1, function(x) {
            x[x <= 0] <- min(x[x > 0])
            x
        }))
        pathway_metabolic_obj <- SetAssayData(pathway_metabolic_obj, new.data = assay_data)
        pathway_number_weight <- 1 / gene_pathway_number[keep, ]
        #
        mean_exp_eachCellType <- t(AverageExpression(pathway_metabolic_obj)$RNA)
        ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
        # exclude the extreme ratios
        col_quantile <- apply(ratio_exp_eachCellType, 2, function(x) quantile(x, na.rm = T))
        col_q1 <- col_quantile["25%", ]
        col_q3 <- col_quantile["75%", ]
        col_upper <- col_q3 * 3
        col_lower <- col_q1 / 3
        outliers <- apply(ratio_exp_eachCellType, 1, function(x) {
            any((x > col_upper) | (x < col_lower))
        })

        if (sum(!outliers) < 3) next

        keep <- names(outliers)[!outliers]
        pathway_metabolic_obj <- subset(pathway_metabolic_obj, features = keep)
        pathway_number_weight <- 1 / gene_pathway_number[keep, ]
        mean_exp_eachCellType <- t(AverageExpression(pathway_metabolic_obj)$RNA)
        ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
        mean_exp_pathway <- apply(ratio_exp_eachCellType, 2, function(x) weighted.mean(x, pathway_number_weight / sum(pathway_number_weight)))
        mean_expression_shuffle[p, ] <- mean_exp_pathway[cell_types]
        mean_expression_noshuffle[p, ] <- mean_exp_pathway[cell_types]
        pathway_metabolic_data <- GetAssayData(pathway_metabolic_obj)

        ## shuffle 5000 times:
        ## define the functions
        group_mean <- function(x) {
            sapply(
                cell_types,
                function(y) rowMeans(pathway_metabolic_data[, shuffle_cell_types_list[[x]] == y, drop = F])
            )
        }
        column_weigth_mean <- function(x) {
            apply(ratio_exp_eachCellType_list[[x]], 2, function(y) weighted.mean(y, weight_values))
        }
        #####
        times <- 1:ntimes
        weight_values <- pathway_number_weight / sum(pathway_number_weight)
        # shuffle_cell_types_list <- mclapply(times, function(x) sample(all_cell_types), mc.cores = ncores)
        shuffle_cell_types_list <- lapply(times, function(x) sample(all_cell_types))
        names(shuffle_cell_types_list) <- times
        # mean_exp_eachCellType_list <- mclapply(times, function(x) group_mean(x), mc.cores = ncores)
        mean_exp_eachCellType_list <- lapply(times, function(x) group_mean(x))
        # ratio_exp_eachCellType_list <- mclapply(times, function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]), mc.cores = ncores)
        ratio_exp_eachCellType_list <- lapply(times, function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
        # mean_exp_pathway_list <- mclapply(times, function(x) column_weigth_mean(x), mc.cores = ncores)
        mean_exp_pathway_list <- lapply(times, function(x) column_weigth_mean(x))

        shuffle_results <- matrix(unlist(mean_exp_pathway_list), ncol = length(cell_types), byrow = T)
        rownames(shuffle_results) <- times
        colnames(shuffle_results) <- cell_types
        for (c in cell_types) {
            if (is.na(mean_expression_shuffle[p, c])) next
            if (mean_expression_shuffle[p, c] > 1) {
                pval <- sum(shuffle_results[, c] > mean_expression_shuffle[p, c]) / ntimes
            } else if (mean_expression_shuffle[p, c] < 1) {
                pval <- sum(shuffle_results[, c] < mean_expression_shuffle[p, c]) / ntimes
            }
            if (pval > 0.01) mean_expression_shuffle[p, c] <- NA ### NA is  blank in heatmap
            pvalues_mat[p, c] <- pval
        }
    }
    all_NA <- rowAlls(is.na(mean_expression_shuffle))
    mean_expression_shuffle <- mean_expression_shuffle[!all_NA, , drop = F]
    # heatmap
    dat <- mean_expression_shuffle

    sort_row <- c()
    sort_column <- c()

    for (i in colnames(dat)) {
        select_row <- which(rowMaxs(dat, na.rm = T) == dat[, i])
        tmp <- rownames(dat)[select_row][order(dat[select_row, i], decreasing = T)]
        sort_row <- c(sort_row, tmp)
    }
    sort_column <- apply(dat[sort_row, , drop = F], 2, function(x) order(x)[nrow(dat)])
    sort_column <- names(sort_column)
    dat[is.na(dat)] <- 1

    heatmapfile <- file.path(subset_dir, "KEGGpathway_activity_heatmap.png")
    hmdata <- dat[sort_row, sort_column, drop = F]
    cnames <- sapply(colnames(hmdata), function(x) {paste0(grouping_prefix, x)})
    colnames(hmdata) <- cnames
    hm_devpars = heatmap_devpars
    if (is.null(hm_devpars$res)) {
        hm_devpars$res = 100
    }
    if (is.null(hm_devpars$width)) {
        hm_devpars$width = 300 + max(nchar(rownames(hmdata))) * 8 + ncol(hmdata) * 12
    }
    if (is.null(hm_devpars$height)) {
        hm_devpars$height = 400 + max(nchar(colnames(hmdata))) * 8 + nrow(hmdata) * 12
    }
    plotHeatmap(
        hmdata,
        args = list(
            name = "Pathway activity",
            rect_gp = gpar(col = "white", lwd = 0.5),
            row_names_side = "left",
            row_dend_side = "right",
            row_names_max_width = max_text_width(
                rownames(hmdata),
                gp = gpar(fontsize = 12)
            ),
            row_dend_width = unit(30, "mm"),
            cluster_columns = FALSE
        ),
        devpars = hm_devpars,
        outfile = heatmapfile
    )


    write.table(mean_expression_noshuffle, file = file.path(subset_dir, "KEGGpathway_activity_noshuffle.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
    write.table(mean_expression_shuffle, file = file.path(subset_dir, "KEGGpathway_activity_shuffle.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
    write.table(pvalues_mat, file = file.path(subset_dir, "KEGGpathway_activity_shuffle_pvalue.txt"), row.names = T, col.names = T, quote = F, sep = "\t")

    # boxplot show the distribution of pathway activity
    scRNA_dat <- as.data.frame(mean_expression_noshuffle)
    scRNA_dat$X <- NULL

    scRNA_df <- melt(as.matrix(scRNA_dat))
    scRNA_df <- scRNA_df[!is.na(scRNA_df$value), ]
    colnames(scRNA_df)[ncol(scRNA_df) - 1] <- "variable"
    scRNA_df$variable <- sapply(scRNA_df$variable, function(x) {paste0(grouping_prefix, x)})
    violinfile <- file.path(subset_dir, "pathway_activity_violinplot.png")
    vio_devpars = violin_devpars
    if (is.null(vio_devpars$res)) {
        vio_devpars$res = 100
    }
    if (is.null(vio_devpars$width)) {
        vio_devpars$width = 100 + ncol(scRNA_df) * 50
    }
    if (is.null(hm_devpars$height)) {
        vio_devpars$height = 1000
    }
    plotViolin(
        scRNA_df,
        args = list(
            mapping = aes(x = variable, y = value, fill = variable),
            trim = F,
            size = 0.2,
            show.legend = F,
            width = 1.2
        ),
        ggs = c(
            "scale_y_continuous(limits = c(0, 3), breaks = 0:3, labels = 0:3)",
            'labs(y = "Metabolic Pathway Activity", x=NULL)',
            'stat_summary(
                aes(x = variable, y = value),
                fun = median,
                geom = "point",
                size = 1,
                color = "black"
            )',
            "theme_prism(axis_text_angle = 90)"
        ),
        devpars = vio_devpars,
        outfile = violinfile
    )
}

if (is.null(subsetting)) {
    do_one_subset(NULL)
} else {
    subsets <- unique(sobj@meta.data[[subsetting]])

    if (ncores == 1) {
        lapply(subsets, do_one_subset)
    } else {
        x <- mclapply(subsets, do_one_subset, mc.cores = ncores)
        if (any(unlist(lapply(x, class)) == "try-error")) {
            stop("mclapply error")
        }
    }
}
