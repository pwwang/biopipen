library(rlang)
library(parallel)
library(matrixStats)
library(enrichit)
library(Seurat)
library(biopipen.utils)
library(plotthis)
library(tidyseurat)

sobjfile <- {{ in.sobjfile | r }}
outdir <- {{ out.outdir | r }}
ntimes <- {{ envs.ntimes | r }}
ncores <- {{ envs.ncores | r }}
gmtfile <- {{ envs.gmtfile | r }}
subset_by <- {{ envs.subset_by | r }}
group_by <- {{ envs.group_by | r }}
plots <- {{ envs.plots | r }}
cases <- {{ envs.cases | r }}

set.seed(8525)

log <- get_logger()
reporter <- get_reporter()

log$info("Loading Seurat object ...")
sobj <- read_obj(sobjfile)
assay <- DefaultAssay(sobj)

defaults <- list(
    ntimes = ntimes,
    subset_by = subset_by,
    group_by = group_by,
    plots = plots
)
log$info("Expanding cases ...")
default_case <- subset_by %||% "DEFAULT"
cases <- expand_cases(
    cases,
    defaults,
    function(name, case) {
        if (is.null(case$group_by)) {
            stop("'group_by' is required in case: ", name)
        }
        stats::setNames(list(case), name)
    },
    default_case = default_case)

log$info("Loading metabolic pathways ...")
pathways <- ParseGMT(gmtfile)
pathway_names <- names(pathways)
metabolics <- unique(as.vector(unname(unlist(pathways))))

## calculate how many pathways of one gene involved.
num_of_pathways <- function(overlapgenes) {
    filter_pathways <- list()
    for (p in pathway_names) {
        genes <- pathways[[p]]
        common_genes <- intersect(genes, overlapgenes)
        if (length(common_genes >= 5)) {
            filter_pathways[[p]] <- common_genes
        }
    }

    all_genes <- unique(as.vector(unlist(filter_pathways)))
    gene_times <- data.frame(
        num = rep(0, length(all_genes)),
        row.names = all_genes
    )
    for (p in pathway_names) {
        for (g in filter_pathways[[p]]) {
            gene_times[g, "num"] <- gene_times[g, "num"] + 1
        }
    }
    gene_times
}

do_subset <- function(
    object,
    caseinfo,
    subset_by,
    subset_val,
    ntimes,
    group_by,
    plots
) {
    if (!is.null(subset_val)) {
        log$info("- Handling subset: {subset_by} = {subset_val} ...")
        object <- tryCatch(
            filter(object, !!sym(subset_by) == subset_val & !is.na(!!sym(group_by))),
            error = function(e) NULL
        )

        if (is.null(object) || ncol(object) < 5) {
            msg <- paste0("  ! skipped. Subset has less than 5 cells: ", subset_by, " = ", subset_val)
            log$warn(msg)
            reporter$add(list(kind = "error", content = msg), h1 = caseinfo$name)
            return(NULL)
        }
    }

    all_groups <- object@meta.data[[group_by]]
    if (!is.factor(all_groups)) {
        all_groups <- factor(all_groups)
    }
    # order by levels(all_groups)
    groups <- intersect(levels(all_groups), unique(all_groups))

    gene_pathway_number <- num_of_pathways(intersect(rownames(object), metabolics))

    ## Calculate the pathway activities
    # mean ratio of genes in each pathway for each cell type
    mean_expression_shuffle <- matrix(
        NA,
        nrow = length(pathway_names),
        ncol = length(groups),
        dimnames = list(pathway_names, groups)
    )
    mean_expression_noshuffle <- matrix(
        NA,
        nrow = length(pathway_names),
        ncol = length(groups),
        dimnames = list(pathway_names, groups)
    )
    ### calculate the pvalues using shuffle method
    pvalues_mat <- matrix(
        NA,
        nrow = length(pathway_names),
        ncol = length(groups),
        dimnames = (list(pathway_names, groups))
    )

    for (pi in seq_along(pathway_names)) {
        p <- pathway_names[pi]
        log$info("  Pathway ({pi}/{length(pathway_names)}): {p} ...")
        genes <- pathways[[p]]
        genes_comm <- intersect(genes, rownames(object))
        # genes_expressed <- names(rowSums(object)[rowSums(object) > 0])
        # genes_comm <- intersect(genes_comm, genes_expressed)
        if (length(genes_comm) < 5) next

        # Errored if default assay is SCT
        # Issue with Seurat?
        # pathway_metabolic_obj <- subset(object, features = genes_comm)
        # assay <- DefaultAssay(object)
        ## AggregateExpression raises Warning: The counts layer for the integrated assay is empty. Skipping assay.
        mean_exp_eachCellType <- suppressMessages(AverageExpression(object, features = genes_comm, assays = assay, group.by = group_by))[[assay]]

        # remove genes which are zeros in any celltype to avoid extreme ratio value
        keep <- rownames(mean_exp_eachCellType)[rowAlls(as.matrix(mean_exp_eachCellType) > 0.001, useNames = F)]
        if (length(keep) < 3) next

        # using the loweset value to replace zeros for avoiding extreme ratio value
        # pathway_metabolic_obj <- subset(object, features = keep)
        assay_data = GetAssayData(object, assay = assay, layer = "data")[keep, , drop = F]
        assay_data <- t(apply(assay_data, 1, function(x) {
            x[x <= 0] <- min(x[x > 0])
            x
        }))
        pathway_metabolic_obj <- suppressWarnings(CreateSeuratObject(CreateAssayObject(data = assay_data), assay = assay))
        pathway_metabolic_obj[[group_by]] <- object[[group_by]]
        Idents(pathway_metabolic_obj) <- Idents(object)
        pathway_number_weight <- 1 / gene_pathway_number[keep, ]
        #
        mean_exp_eachCellType <- t(suppressMessages(AverageExpression(pathway_metabolic_obj, assays = assay, group.by = group_by)[[assay]]))
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
        pathway_metabolic_obj <- suppressWarnings(subset(pathway_metabolic_obj, features = keep))
        pathway_number_weight <- 1 / gene_pathway_number[keep, ]
        mean_exp_eachCellType <- t(suppressMessages(AverageExpression(pathway_metabolic_obj, assays = assay, group.by = group_by)[[assay]]))
        ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
        mean_exp_pathway <- apply(ratio_exp_eachCellType, 2, function(x) weighted.mean(x, pathway_number_weight / sum(pathway_number_weight)))
        mean_expression_shuffle[p, ] <- mean_exp_pathway[groups]
        mean_expression_noshuffle[p, ] <- mean_exp_pathway[groups]
        pathway_metabolic_data <- GetAssayData(pathway_metabolic_obj)

        ## shuffle 5000 times:
        ## define the functions
        group_mean <- function(x) {
            sapply(
                groups,
                function(y) rowMeans(pathway_metabolic_data[, shuffle_groups_list[[x]] == y, drop = F])
            )
        }
        column_weigth_mean <- function(x) {
            apply(ratio_exp_eachCellType_list[[x]], 2, function(y) weighted.mean(y, weight_values))
        }
        #####
        times <- 1:ntimes
        weight_values <- pathway_number_weight / sum(pathway_number_weight)
        shuffle_groups_list <- mclapply(times, function(x) sample(all_groups), mc.cores = ncores)
        # shuffle_groups_list <- lapply(times, function(x) sample(all_groups))
        names(shuffle_groups_list) <- times
        mean_exp_eachCellType_list <- mclapply(times, function(x) group_mean(x), mc.cores = ncores)
        # mean_exp_eachCellType_list <- lapply(times, function(x) group_mean(x))
        ratio_exp_eachCellType_list <- mclapply(times, function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]), mc.cores = ncores)
        # ratio_exp_eachCellType_list <- lapply(times, function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
        mean_exp_pathway_list <- mclapply(times, function(x) column_weigth_mean(x), mc.cores = ncores)
        # mean_exp_pathway_list <- lapply(times, function(x) column_weigth_mean(x))

        shuffle_results <- matrix(unlist(mean_exp_pathway_list), ncol = length(groups), byrow = T)
        rownames(shuffle_results) <- times
        colnames(shuffle_results) <- groups
        for (c in groups) {
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
    all_NA <- rowAlls(is.na(as.matrix(mean_expression_shuffle)), useNames = F)
    if (all(all_NA)) {
        log$warn("  ! All pathways are NA after shuffling.")
        # keep at least 3 pathways for plotting
        mean_expression_shuffle <- mean_expression_shuffle[1:3, , drop = F]
        mean_expression_shuffle[is.na(mean_expression_shuffle)] <- 1
    } else {
        mean_expression_shuffle <- mean_expression_shuffle[!all_NA, , drop = F]
    }
    # heatmap
    dat <- mean_expression_shuffle

    sort_row <- c()
    sort_column <- c()

    for (i in colnames(dat)) {
        select_row <- which(rowMaxs(dat, na.rm = TRUE, useNames = FALSE) == dat[, i])
        tmp <- rownames(dat)[select_row][order(dat[select_row, i], decreasing = TRUE)]
        sort_row <- unique(c(sort_row, tmp))
    }
    sort_column <- apply(dat[sort_row, , drop = FALSE], 2, function(x) order(x)[nrow(dat)])
    sort_column <- names(sort_column)
    dat[is.na(dat)] <- 1
    dat <- dat[sort_row, sort_column, drop = FALSE]

    if (!is.null(subset_by)) {
        prefix <- file.path(caseinfo$prefix, paste0(slugify(subset_by), "_", slugify(subset_val), "."))
        h2 <- paste0(subset_by, ": ", subset_val)
    } else if (length(cases) > 1) {
        prefix <- paste0(caseinfo$prefix, "/No_Subsetting/")
        dir.create(prefix, showWarnings = FALSE, recursive = TRUE)
        h2 <- "No Subsetting"
    } else {
        prefix <- paste0(caseinfo$prefix, "/")
        h2 <- "#"
    }

    write.table(
        mean_expression_noshuffle,
        file = paste0(prefix, "pathway_activity_noshuffle.txt"),
        row.names = TRUE,
        col.names = TRUE,
        quote = FALSE,
        sep = "\t"
    )
    write.table(
        mean_expression_shuffle,
        file = paste0(prefix, "pathway_activity_shuffle.txt"),
        row.names = TRUE,
        col.names = TRUE,
        quote = FALSE,
        sep = "\t"
    )
    write.table(pvalues_mat,
        file = paste0(prefix, "pathway_activity_shuffle_pvalue.txt"),
        row.names = TRUE,
        col.names = TRUE,
        quote = FALSE,
        sep = "\t"
    )

    for (plotname in names(plots)) {
        plotargs <- plots[[plotname]]
        plotargs$devpars <- plotargs$devpars %||% list()
        plotargs <- extract_vars(plotargs, "devpars", "plot_type")
        devpars <- devpars %||% list()
        devpars$res <- devpars$res %||% 100
        if (plot_type == "merged_heatmap") { next }
        log$info("  Plotting: {plotname} ...")
        if (plot_type %in% c("violin", "box", "boxplot")) {
            plotfn <- if (plot_type == "violin") plotthis::ViolinPlot else plotthis::BoxPlot
            # boxplot show the distribution of pathway activity
            scRNA_dat <- as.data.frame(mean_expression_noshuffle)
            scRNA_dat$X <- NULL

            # scRNA_df <- reshape2::melt(as.matrix(scRNA_dat))
            # scRNA_df <- scRNA_df[!is.na(scRNA_df$value), ]
            # colnames(scRNA_df)[ncol(scRNA_df) - 1] <- "variable"
            scRNA_dat$Pathways <- rownames(scRNA_dat)
            scRNA_dat <- tidyr::pivot_longer(
                scRNA_dat,
                cols = -c(Pathways),
                names_to = group_by,
                values_to = "Pathway Activity"
            )

            plotargs$data <- scRNA_dat
            plotargs$x <- group_by
            plotargs$y <- "Pathway Activity"
            plotargs$keep_empty <- TRUE

            p <- do_call(plotfn, plotargs)
            devpars$width <- devpars$width %||% (attr(p, "width") * 2 * devpars$res) %||% 1000
            devpars$height <- devpars$height %||% (attr(p, "height") * 2 * devpars$res) %||% 1000
        } else {  # heatmap
            minval <- min(dat)
            maxval <- max(dat)
            dis <- max(1 - minval, maxval - 1)
            minval <- 1 - dis
            maxval <- 1 + dis
            dat <- as.data.frame(t(dat))  # rows: groups, columns: pathways
            dat[[group_by]] <- rownames(dat)
            plotargs$data <- dat
            plotargs$columns_by <- group_by
            plotargs$in_form <- "wide-rows"
            plotargs$name <- plotargs$name %||% "Pathway Activity"
            plotargs$rows_name <- plotargs$rows_name %||% "Pathways"
            plotargs$show_row_names <- plotargs$show_row_names %||% TRUE
            plotargs$lower_cutoff <- plotargs$lower_cutoff %||% minval
            plotargs$upper_cutoff <- plotargs$upper_cutoff %||% maxval
            plotargs$row_name_annotation <- plotargs$row_name_annotation %||% FALSE
            plotargs$row_names_side <- plotargs$row_names_side %||% "left"
            plotargs$show_column_names <- plotargs$show_column_names %||% TRUE

            p <- do_call(plotthis::Heatmap, plotargs)
            devpars$width <- devpars$width %||% (attr(p, "width") * devpars$res) %||% 1000
            devpars$height <- devpars$height %||% (attr(p, "height") * devpars$res) %||% 1000
        }

        plotprefix <- paste0(prefix, slugify(plotname))
        png(paste0(plotprefix, ".png"), res = devpars$res, width = devpars$width, height = devpars$height)
        print(p)
        dev.off()

        descr <- plotargs$descr %||% paste0(
            plotname, " a ", plotargs$plot_type, " plot of pathway activity for ", group_by, ". "
        )

        reporter$add(
            list(name = plotname, contents = list(
                list(kind = "descr", content = descr),
                reporter$image(plotprefix, c(), FALSE))
            ),
            h1 = caseinfo$name,
            h2 = h2,
            ui = "tabs"
        )
    }

    return(dat)
}


do_case <- function(casename) {
    log$info("Processing case: {casename} ...")
    case <- cases[[casename]]
    if (is.null(case) || length(case) == 0) {
        log$warn("  Case skipped.")
        return(NULL)
    }
    caseinfo <- case_info(casename, outdir, create = TRUE)

    if (is.null(case$subset_by)) {
        result <- do_subset(
            sobj,
            caseinfo = caseinfo,
            subset_by = NULL,
            subset_val = NULL,
            ntimes = case$ntimes,
            group_by = case$group_by,
            plots = case$plots
        )
    } else {
        sobj_avail <- filter(sobj, !is.na(!!sym(case$subset_by)))
        if (ncol(sobj_avail) < 5) {
            stop("Not enough cells (< 5) for subset: ", case$subset_by)
        }

        subsets <- unique(sobj@meta.data[[case$subset_by]])
        result <- NULL
        for (ss in subsets) {
            tmp <- do_subset(
                sobj_avail,
                caseinfo = caseinfo,
                subset_by = case$subset_by,
                subset_val = ss,
                ntimes = case$ntimes,
                group_by = case$group_by,
                plots = case$plots
            )
            if (is.null(tmp)) { next }
            tmp[[case$group_by]] <- rownames(tmp)
            tmp[[case$subset_by]] <- ss
            rownames(tmp) <- NULL
            if (is.null(result)) {
                result <- tmp
            } else {
                all_columns <- union(colnames(result), colnames(tmp))
                result[, setdiff(all_columns, colnames(result))] <- 1
                tmp[, setdiff(all_columns, colnames(tmp))] <- 1
                result <- rbind(result, tmp)
            }
        }
        uniq_subsets <- unique(result[[case$subset_by]])
        result[[case$subset_by]] <- factor(
            result[[case$subset_by]],
            levels = if (is.factor(sobj@meta.data[[case$subset_by]])) {
                intersect(levels(sobj@meta.data[[case$subset_by]]), uniq_subsets)
            } else {
                uniq_subsets
            }
        )
    }
    uniq_groups <- unique(result[[case$group_by]])
    result[[case$group_by]] <- factor(
        result[[case$group_by]],
        levels = if (is.factor(sobj@meta.data[[case$group_by]])) {
            intersect(levels(sobj@meta.data[[case$group_by]]), uniq_groups)
        } else {
            uniq_groups
        }
    )

    for (plotname in names(case$plots)) {
        plotargs <- case$plots[[plotname]]
        if (is.null(plotargs$plot_type)) {
            stop("'plot_type' is required in plot args: ", plotname, " in case: ", casename)
        }
        plotargs$devpars <- plotargs$devpars %||% list()
        plotargs <- extract_vars(plotargs, "devpars", "plot_type")
        if (plot_type != "merged_heatmap") {
            next
        }
        log$info("  Plotting: {plotname} ...")

        plotargs$data <- result
        plotargs$name <- plotargs$name %||% "Pathway Activity"
        plotargs$in_form <- "wide-rows"
        plotargs$columns_by <- case$group_by
        plotargs$show_row_names <- plotargs$show_row_names %||% TRUE
        minval <- min(as.matrix(result[, setdiff(colnames(result), c(case$group_by, case$subset_by))]))
        maxval <- max(as.matrix(result[, setdiff(colnames(result), c(case$group_by, case$subset_by))]))
        dis <- max(1 - minval, maxval - 1)
        minval <- 1 - dis
        maxval <- 1 + dis
        plotargs$lower_cutoff <- plotargs$lower_cutoff %||% minval
        plotargs$upper_cutoff <- plotargs$upper_cutoff %||% maxval
        plotargs$row_name_annotation <- plotargs$row_name_annotation %||% FALSE
        plotargs$row_names_side <- plotargs$row_names_side %||% "left"
        plotargs$show_column_names <- plotargs$show_column_names %||% TRUE

        if (!is.null(case$subset_by)) {
            plotargs$columns_split_by <- case$subset_by
        }
        p <- do_call(plotthis::Heatmap, plotargs)

        devpars <- devpars %||% list()
        devpars$res <- devpars$res %||% 100
        devpars$width <- devpars$width %||% (attr(p, "width") * devpars$res) %||% 1000
        devpars$height <- devpars$height %||% (attr(p, "height") * devpars$res) %||% 1000

        prefix <- file.path(caseinfo$prefix, paste0(slugify(plotname), ".merged_heatmap"))
        png(paste0(prefix, ".png"), res = devpars$res, width = devpars$width, height = devpars$height)
        print(p)
        dev.off()

        descr <- plotargs$descr %||% "Merged Heatmaps for Pathway Activity of all subsets."

        reporter$add(
            list(name = plotname, contents = list(
                list(kind = "descr", content = descr),
                reporter$image(prefix, c(), FALSE)
            )),
            h1 = casename,
            h2 = "Merged Heatmaps",
            ui = "tabs"
        )
    }
}

sapply(names(cases), do_case)

reporter$save(dirname(outdir))
