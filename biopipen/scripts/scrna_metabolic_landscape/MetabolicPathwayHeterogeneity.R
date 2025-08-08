library(gtools)
library(rlang)
library(Matrix)
library(sparseMatrixStats)
library(Seurat)
library(tidyseurat)
library(biopipen.utils)

sobjfile <- {{ in.sobjfile | r }}
outdir <- {{ out.outdir | r }}
joboutdir <- {{ job.outdir | r }}
gmtfile <- {{ envs.gmtfile | r }}
select_pcs <- {{ envs.select_pcs | r }}
ncores <- {{ envs.ncores | r }}
pathway_pval_cutoff <- {{ envs.pathway_pval_cutoff | r }}
subset_by <- {{ envs.subset_by | r }}
group_by <- {{ envs.group_by | r }}
fgsea_args <- {{ envs.fgsea_args | r }}
plots <- {{ envs.plots | r }}
cases <- {{ envs.cases | r }}

set.seed(8525)

log <- get_logger()
reporter <- get_reporter()

log$info("Loading Seurat object ...")
sobj <- read_obj(sobjfile)

defaults <- list(
    subset_by = subset_by,
    group_by = group_by,
    fgsea_args = fgsea_args,
    plots = plots,
    select_pcs = select_pcs,
    pathway_pval_cutoff = pathway_pval_cutoff
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


do_subset <- function(object, caseinfo, subset_by, subset_val, group_by, plots, select_pcs, pathway_pval_cutoff) {
    if (!is.null(subset_by)) {
        log$info("- Handling subset: {subset_by} = {subset_val} ...")
        object <- tryCatch(
            filter(object, !!sym(subset_by) == subset_val & !is.na(!!sym(group_by))),
            error = function(e) NULL
        )
    }
    if (!is.null(subset_by)) {
        h1 <- paste0(subset_by, ": ", subset_val)
        h2 <- group_by
        odir <- file.path(caseinfo$prefix, slugify(paste0(subset_by, "_", subset_val)))
    } else if (length(cases) > 1) {
        h1 <- "No Subsetting"
        h2 <- group_by
        odir <- file.path(caseinfo$prefix, "No_Subsetting")
    } else {
        h1 <- group_by
        h2 <- "#"
        odir <- caseinfo$prefix
    }
    if (is.null(object) || ncol(object) < 5) {
        msg <- paste0("  ! skipped. Subset has less than 5 cells: ", subset_by, " = ", subset_val)
        log$warn(msg)
        reporter$add(list(kind = "error", content = msg), h1 = h1, h2 = h2)
        return(NULL)
    }

    dir.create(odir, showWarnings = FALSE)

    features <- intersect(rownames(object), metabolics)
    groups <- unique(as.character(object@meta.data[[group_by]]))

    enrich_data_df <- NULL
    pc_plotdata <- data.frame(
        x = numeric(),
        y = numeric(),
        sel = character(),
        group = character()
    )

    for (group in groups) {
        log$info("  {group_by}: {group} ...")
        each_metabolic_obj <- subset(object, subset = !!sym(group_by) == group)
        if (ncol(each_metabolic_obj) < 5) {
            log$warn("  ! skipped. Group has less than 5 cells: {group}")
            next()
        }
        each_metabolic_exprs <- GetAssayData(each_metabolic_obj)[features, , drop = FALSE]
        each_metabolic_exprs <- each_metabolic_exprs[rowSums(each_metabolic_exprs) > 0, , drop=FALSE]
        if (ncol(each_metabolic_obj) < 5) {
            log$warn("  ! skipped. Group has less than 5 active cells: {group}")
            next()
        }
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
        pre_rank_matrix <- as.matrix(rowSums(abs(pca$rotation[, 1:selected_pcs, drop = FALSE])))
        pre_rank_matrix <- unlist(as.list(as.data.frame(t(pre_rank_matrix))))

        fgsea_args <- fgsea_args %||% list()
        fgsea_args$ranks <- pre_rank_matrix
        fgsea_args$genesets <- pathways
        fgsea_args$nproc <- fgsea_args$nproc %||% ncores

        tmp <- do_call(RunGSEA, fgsea_args)
        tmp[[group_by]] <- group

        if (is.null(enrich_data_df)) {
            enrich_data_df <- tmp
        } else {
            enrich_data_df <- rbind(enrich_data_df, tmp)
        }
    }

    # remove pvalue < 0.01 pathways
    min_pval <- by(enrich_data_df$pval, enrich_data_df$pathway, FUN = min)
    select_pathways <- names(min_pval)[(min_pval <= pathway_pval_cutoff)]
    select_enrich_data_df <- enrich_data_df[enrich_data_df$pathway %in% select_pathways, ]
    # converto pvalue to -log10
    pvals <- select_enrich_data_df$pval
    pvals[pvals <= 0] <- 1e-10
    select_enrich_data_df$pval <- -log10(pvals)

    # sort
    pathway_pv_sum <- by(select_enrich_data_df$pval, select_enrich_data_df$pathway, FUN = sum)
    pathway_order <- names(pathway_pv_sum)[order(pathway_pv_sum, decreasing = T)]
    ########################### top 10
    pathway_order <- pathway_order[1:10]
    select_enrich_data_df <- select_enrich_data_df[select_enrich_data_df$pathway %in% pathway_order, ]
    ########################################
    select_enrich_data_df[[group_by]] <- factor(select_enrich_data_df[[group_by]], levels = gtools::mixedsort(groups))
    select_enrich_data_df$pathway <- factor(select_enrich_data_df$pathway, levels = pathway_order)

    write.table(
        as.data.frame(select_enrich_data_df),
        file = file.path(odir, "pathway_heterogeneity.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    for (plot in names(plots)) {
        plotargs <- plots[[plot]]
        plotargs$devpars <- plotargs$devpars %||% list()
        plotargs$devpars$res <- plotargs$devpars$res %||% 100

        if (plotargs$plot_type == "dot") {
            plotargs$x <- plotargs$x %||% group_by
            plotargs$y <- plotargs$y %||% "pathway"
            plotargs$fill_by <- plotargs$fill_by %||% "NES"
            plotargs$size_by <- plotargs$size_by %||% "pval"
            plotargs$add_bg <- plotargs$add_bg %||% TRUE
            plotargs$x_text_angle <- plotargs$x_text_angle %||% 90
            plotfn <- plotthis::DotPlot
        } else {
            stop("Unknown plot type: ", plotargs$plot_type)
        }

        p <- do_call(plotfn, c(list(select_enrich_data_df), plotargs))
        plotprefix <- file.path(odir, slugify(plot))
        plotargs$devpars$width <- plotargs$devpars$width %||% (attr(p, "width") * plotargs$devpars$res) %||% 800
        plotargs$devpars$height <- plotargs$devpars$height %||% (attr(p, "height") * plotargs$devpars$res) %||% 600
        plotargs$devpars$height <- max(plotargs$devpars$height, plotargs$devpars$width / 1.5)
        png(
            filename = paste0(plotprefix, ".png"),
            width = plotargs$devpars$width,
            height = plotargs$devpars$height,
            res = plotargs$devpars$res
        )
        print(p)
        dev.off()

        reporter$add(
            list(
                name = plot,
                contents = list(
                    list(kind = "descr", content = plotargs$descr %||% plot),
                    reporter$image(plotprefix, c(), FALSE, kind = "image")
                )
            ),
            h1 = h1,
            h2 = h2,
            ui = "tabs"
        )
    }
}


do_case <- function(casename) {
    log$info("Processing case: {casename} ...")
    case <- cases[[casename]]
    caseinfo <- case_info(casename, outdir, create = TRUE)

    if (is.null(case$subset_by)) {
        result <- do_subset(
            sobj,
            caseinfo = caseinfo,
            subset_by = NULL,
            subset_val = NULL,
            group_by = case$group_by,
            plots = case$plots,
            select_pcs = case$select_pcs,
            pathway_pval_cutoff = case$pathway_pval_cutoff
        )
    } else {
        sobj_avail <- filter(sobj, !is.na(!!sym(case$subset_by)))
        subsets <- unique(sobj@meta.data[[case$subset_by]])

        lapply(
            subsets,
            function(ss) {
                do_subset(
                    sobj_avail,
                    caseinfo = caseinfo,
                    subset_by = case$subset_by,
                    subset_val = ss,
                    group_by = case$group_by,
                    plots = case$plots,
                    select_pcs = case$select_pcs,
                    pathway_pval_cutoff = case$pathway_pval_cutoff
                )
            }
        )
    }
}

sapply(names(cases), do_case)

reporter$save(dirname(outdir))
