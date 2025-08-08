library(rlang)
library(Seurat)
library(biopipen.utils)
library(enrichit)
library(tidyseurat)

sobjfile <- {{ in.sobjfile | r }}
outdir <- {{ out.outdir | r }}
ncores <- {{ envs.ncores | r }}
prerank_method <- {{ envs.prerank_method | r }}
gmtfile <- {{ envs.gmtfile | r }}
subset_by <- {{ envs.subset_by | r }}
group_by <- {{ envs.group_by | r }}
comparisons <- {{ envs.comparisons | r }}
fgsea_args <- {{ envs.fgsea_args | r }}
plots <- {{ envs.plots | r }}
cases <- {{ envs.cases | r }}

set.seed(8525)

log <- get_logger()
reporter <- get_reporter()

log$info("Loading Seurat object ...")
sobj <- read_obj(sobjfile)

defaults <- list(
    prerank_method = prerank_method,
    subset_by = subset_by,
    group_by = group_by,
    comparisons = comparisons,
    fgsea_args = fgsea_args,
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


do_comparison <- function(object, caseinfo, subset_by, subset_val, group_by, group1, group2, prerank_method, plots, fgsea_args) {
    log$info("  {group_by}: {group1} vs {group2} ...")
    if (!is.null(group2)) {
        # object <- subset(object, subset = !!sym(group_by) %in% c(group1, group2))
        object <- tryCatch(
            filter(object, !!sym(group_by) %in% c(group1, group2)),
            error = function(e) NULL
        )
    }

    if (!is.null(subset_by)) {
        if (length(cases) == 1 && identical(caseinfo$name, subset_by)) {
            # No need to show the case name in report
            h1 <- paste0(subset_by, ": ", subset_val)
            h2 <- paste0(group_by, ": ", group1, " vs ", group2 %||% "REST")
            h3 <- "#"
        } else {
            h1 <- caseinfo$name
            h2 <- paste0(subset_by, ": ", subset_val)
            h3 <- paste0(group_by, ": ", group1, " vs ", group2 %||% "REST")
        }
        odir <- file.path(caseinfo$prefix, slugify(paste0(subset_by, "_", subset_val)))
    } else if (length(cases) > 1) {
        h1 <- caseinfo$name
        h2 <- "#"
        h3 <- paste0(group_by, ": ", group1, " vs ", group2 %||% "REST")
        odir <- file.path(caseinfo$prefix, "No_Subsetting")
    } else {
        h1 <- paste0(group_by, ": ", group1, " vs ", group2 %||% "REST")
        h2 <- "#"
        h3 <- "#"
        odir <- caseinfo$prefix
    }

    if (is.null(object) || ncol(object) < 10) {
        msg <- paste0("  ! skipped. Groups together have less than 10 cells: ", group_by, " = ", group1, " vs ", group2)
        log$warn(msg)
        reporter$add(
            list(kind = "error", content = msg),
            h1 = h1,
            h2 = h2,
            h3 = h3
        )
        return(invisible())
    }

    classes <- as.character(object@meta.data[[group_by]])
    if (!group1 %in% classes) {
        stop("Group '", group1, "' not found in '", group_by, "' column of the Seurat object.")
    }
    if (!is.null(group2) && !group2 %in% classes) {
        stop("Group '", group2, "' not found in '", group_by, "' column of the Seurat object.")
    }
    classes[classes != group1] <- "Other"
    if (any(table(classes) < 5)) {
        msg <- paste0(
            "  ! skipped. Group has less than 5 cells: ",
            paste(names(table(classes)[table(classes) < 5]), collapse = ", ")
        )
        log$warn(msg)

        reporter$add(
            list(kind = "error", content = msg),
            h1 = h1,
            h2 = h2,
            h3 = h3
        )
        return(invisible())
    }

    features = intersect(rownames(object), metabolics)
    ranks <- RunGSEAPreRank(
        GetAssayData(object)[features, , drop = FALSE],
        classes = object@meta.data[[group_by]],
        case = group1,
        control = group2,
        method = prerank_method
    )

    fgsea_args <- fgsea_args %||% list()
    fgsea_args$ranks <- ranks
    fgsea_args$genesets <- pathways
    fgsea_args$nproc <- fgsea_args$nproc %||% ncores
    result <- do_call(RunGSEA, fgsea_args)

    if (is.null(group2)) {
        odir <- file.path(odir, slugify(paste0(group_by, "_", group1, "_vs_REST")))
    } else {
        odir <- file.path(odir, slugify(paste0(group_by, "_", group1, "_vs_", group2)))
    }

    dir.create(odir, showWarnings = FALSE, recursive = TRUE)
    write.table(as.data.frame(result), file = file.path(odir, "fgsea_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(data.frame(Gene = names(ranks), Rank = ranks), file = file.path(odir, "fgsea_ranks.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

    reporter$add(
        list(kind = "descr", content = "A summary table of the GSEA results"),
        list(kind = "table", src = file.path(odir, "fgsea_results.txt")),
        h1 = h1,
        h2 = h2,
        h3 = h3
    )

    for (plot in names(plots)) {
        plotargs <- plots[[plot]]
        plotargs$level <- plotargs$level %||% "group"
        if (plotargs$level != "group") { next }
        plotargs$devpars <- plotargs$devpars %||% list()
        plotargs$devpars$res <- plotargs$devpars$res %||% 100

        if (identical(plotargs$plot_type, "summary")) {
            p <- do_call(VizGSEA, c(list(result), plotargs))
            plotprefix <- file.path(odir, slugify(plot))
            plotargs$devpars$width <- plotargs$devpars$width %||% (attr(p, "width") * plotargs$devpars$res) %||% 800
            plotargs$devpars$height <- plotargs$devpars$height %||% (attr(p, "height") * plotargs$devpars$res) %||% 600
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
                        reporter$image(plotprefix, "png", FALSE, kind = "image")
                    )
                ),
                h1 = h1,
                h2 = h2,
                h3 = h3,
                ui = "tabs"
            )
        } else {
            plotargs$combine = FALSE
            plotargs$top_term = plotargs$top_term %||% 10
            plotargs$gs <- result$pathway[1:plotargs$top_term]

            ps <- do_call(VizGSEA, c(list(result), plotargs))
            plotprefix <- file.path(odir, slugify(plot))
            devpars <- plotargs$devpars
            images <- list()
            for (pname in names(ps)) {
                p <- ps[[pname]]
                devpars$width <- devpars$width %||% (attr(p, "width") * devpars$res) %||% 800
                devpars$height <- devpars$height %||% (attr(p, "height") * devpars$res) %||% 600
                prefix <- paste0(plotprefix, ".", slugify(pname))
                images[[length(images) + 1]] <- reporter$image(prefix, c(), FALSE, kind = "table_image")
                png(
                    filename = paste0(prefix, ".png"),
                    width = devpars$width,
                    height = devpars$height,
                    res = devpars$res
                )
                print(p)
                dev.off()
            }

            reporter$add(
                list(
                    name = plot,
                    ui = "table_of_images:2",
                    contents = images
                ),
                h1 = h1,
                h2 = h2,
                h3 = h3,
                ui = "tabs"
            )
        }
    }
    result$comparison <- paste0(group1, " vs ", group2 %||% "REST")
    return(result)
}


do_subset <- function(object, caseinfo, subset_by, subset_val, group_by, comparisons, prerank_method, plots, fgsea_args) {
    if (!is.null(subset_by)) {
        log$info("- Handling subset: {subset_by} = {subset_val} ...")
        # object <- subset(object, subset = !!sym(subset_by) == subset_val)
        object <- tryCatch(
            filter(object, !!sym(subset_by) == subset_val & !is.na(!!sym(group_by))),
            error = function(e) NULL
        )

        if (is.null(object) || ncol(object) < 5) {
            if (length(cases) == 1 && identical(caseinfo$name, subset_by)) {
                # No need to show case name in report
                h1 <- paste0(subset_by, ": ", subset_val)
                h2 <- "#"
            } else {
                h1 <- caseinfo$name
                h2 <- paste0(subset_by, ": ", subset_val)
            }

            msg <- paste0("  ! skipped. Subset has less than 5 cells: ", subset_by, " = ", subset_val)
            log$warn(msg)
            reporter$add(list(kind = "error", content = msg), h1 = h1)
            return(NULL)
        }
    }

    groups <- unique(object@meta.data[[group_by]])
    if (length(comparisons) == 0) {
        result <- do_call(
            rbind, lapply(
                as.character(groups),
                function(group) {
                    do_comparison(object, caseinfo, subset_by, subset_val, group_by, group, NULL, prerank_method, plots, fgsea_args)
                }
            )
        )
    } else {
        result <- do_call(
            rbind, lapply(
                as.character(comparisons),
                function(comparison) {
                    if (grepl(":", comparison)) {
                        group1 <- trimws(unlist(strsplit(comparison, ":")))
                        group2 <- group1[2]
                        group1 <- group1[1]
                    } else {
                        group1 <- comparison
                        group2 <- NULL
                    }
                    do_comparison(object, caseinfo, subset_by, subset_val, group_by, group1, group2, prerank_method, plots, fgsea_args)
                }
            )
        )
    }

    result[["-log10(pval)"]] <- -log10(result$pval)
    result[["-log10(padj)"]] <- -log10(result$padj)

    odir <- NULL
    if (!is.null(subset_by)) {
        if (length(cases) == 1 && identical(caseinfo$name, subset_by)) {
            # No need to show case name in report
            h1 <- paste0(subset_by, ": ", subset_val)
            h2 <- "Summary plots for all comparisons"
            h3 <- "#"
        } else {
            h1 <- caseinfo$name
            h2 <- paste0(subset_by, ": ", subset_val)
            h3 <- "Summary plots for all comparisons"
        }
        odir <- file.path(caseinfo$prefix, slugify(paste0(subset_by, "_", subset_val)))
    } else if (length(cases) > 1) {
        h1 <- caseinfo$name
        h2 <- "Summary plots for all comparisons"
        h3 <- "#"
        odir <- file.path(caseinfo$prefix, "No_Subsetting")
    }

    if (!is.null(odir)) {
        dir.create(odir, showWarnings = FALSE, recursive = TRUE)
    }

    for (plot in names(plots)) {
        plotargs <- plots[[plot]]
        plotargs$level <- plotargs$level %||% "group"
        if (plotargs$level != "subset") { next }
        if (is.null(odir)) {
            stop("'subset_by' is NULL but plot level is 'subset': ", plot, ", use level = 'case' instead.")
        }
        plotargs$devpars <- plotargs$devpars %||% list()
        plotargs$devpars$res <- plotargs$devpars$res %||% 100
        plotargs$plot_type <- plotargs$plot_type %||% "dot"

        if (identical(plotargs$plot_type, "dot")) {
            plotargs$x <- plotargs$x %||% "comparison"
            plotargs$y <- plotargs$y %||% "pathway"
            plotargs$size_by <- plotargs$size_by %||% "NES"
            plotargs$fill_by <- plotargs$fill_by %||% "-log10(padj)"
            plotargs$fill_cutoff <- plotargs$fill_cutoff %||% -log10(0.05)
            plotargs$fill_cutoff_name <- plotargs$fill_cutoff_name %||% "Insignificant"
            plotargs$aspect.ratio <- plotargs$aspect.ratio %||% (length(unique(result$pathway)) / length(unique(result$comparison)) / 4)
            plotargs$x_text_angle <- plotargs$x_text_angle %||% 90
            p <- do_call(plotthis::DotPlot, c(list(result), plotargs))
            plotprefix <- file.path(odir, slugify(plot))
            plotargs$devpars$width <- plotargs$devpars$width %||% (attr(p, "width") * plotargs$devpars$res) %||% 800
            plotargs$devpars$height <- plotargs$devpars$height %||% (attr(p, "height") * plotargs$devpars$res) %||% 600
            png(
                filename = paste0(plotprefix, ".png"),
                width = plotargs$devpars$width,
                height = plotargs$devpars$height,
                res = plotargs$devpars$res
            )
            print(p)
            dev.off()

            reporter$add(
                list(kind = "descr", content = plotargs$descr %||% plot),
                reporter$image(plotprefix, "png", FALSE, kind = "image"),
                h1 = h1,
                h2 = h2,
                h3 = h3
            )
        } else {
            stop("`subset` level plot type not supported yet: ", plotargs$plot_type)
        }
    }
    if (!is.null(subset_by)) {
        result[[subset_by]] <- subset_val
    }

    return(result)
}


do_case <- function(casename) {
    log$info("Processing case: {casename} ...")
    case <- cases[[casename]]
    caseinfo <- case_info(casename, outdir, create = TRUE)

    if (is.null(case$subset_by)) {
        result <- do_subset(
            sobj,
            caseinfo,
            subset_by = NULL,
            subset_val = NULL,
            group_by = case$group_by,
            comparisons = case$comparisons,
            prerank_method = case$prerank_method,
            plots = case$plots,
            fgsea_args = case$fgsea_args
        )
    } else {
        sobj_avail <- filter(sobj, !is.na(!!sym(case$subset_by)))
        subsets <- if (is.factor(sobj_avail@meta.data[[case$subset_by]])) {
            levels(sobj_avail@meta.data[[case$subset_by]])
        } else {
            unique(sobj_avail@meta.data[[case$subset_by]])
        }
        result <- do_call(
            rbind, lapply(
                as.character(subsets),
                function(subset_val) {
                    do_subset(
                        sobj_avail,
                        caseinfo,
                        subset_by = case$subset_by,
                        subset_val = subset_val,
                        group_by = case$group_by,
                        comparisons = case$comparisons,
                        prerank_method = case$prerank_method,
                        plots = case$plots,
                        fgsea_args = case$fgsea_args
                    )
                }
            )
        )
        result[[case$subset_by]] <- factor(result[[case$subset_by]], levels = subsets)
    }
    result$pathway <- factor(result$pathway, levels = pathway_names)

    if (!is.null(case$subset_by)) {
        if (length(cases) == 1 && identical(caseinfo$name, case$subset_by)) {
            h1 <- "Summary plots for all subsets"
            h2 <- "#"
        } else {
            h1 <- caseinfo$name
            h2 <- "Summary plots for all subsets"
        }
    } else if (length(cases) > 1) {
        h1 <- caseinfo$name
        h2 <- "Summary plots for all comparisons"
    } else {
        h1 <- "Summary plots for all comparisons"
        h2 <- "#"
    }

    for (plot in names(plots)) {
        plotargs <- plots[[plot]]
        plotargs$level <- plotargs$level %||% "group"
        if (plotargs$level != "case") { next }
        plotargs$devpars <- plotargs$devpars %||% list()
        plotargs$devpars$res <- plotargs$devpars$res %||% 100
        plotargs$plot_type <- plotargs$plot_type %||% "dot"

        if (identical(plotargs$plot_type, "dot")) {
            plotargs$x <- plotargs$x %||% "comparison"
            plotargs$y <- plotargs$y %||% "pathway"
            plotargs$size_by <- plotargs$size_by %||% "NES"
            plotargs$fill_by <- plotargs$fill_by %||% "-log10(padj)"
            plotargs$fill_cutoff <- plotargs$fill_cutoff %||% -log10(0.05)
            plotargs$fill_cutoff_name <- plotargs$fill_cutoff_name %||% "Insignificant"
            plotargs$aspect.ratio <- plotargs$aspect.ratio %||% (length(unique(result$pathway)) / length(unique(result$comparison)) / 4)
            plotargs$x_text_angle <- plotargs$x_text_angle %||% 90
            if (!is.null(subset_by)) {
                plotargs$split_by <- plotargs$split_by %||% subset_by
            }
            p <- do_call(plotthis::DotPlot, c(list(result), plotargs))
            plotprefix <- file.path(caseinfo$prefix, slugify(plot))
            plotargs$devpars$width <- plotargs$devpars$width %||% (attr(p, "width") * plotargs$devpars$res) %||% 800
            plotargs$devpars$height <- plotargs$devpars$height %||% (attr(p, "height") * plotargs$devpars$res) %||% 600
            png(
                filename = paste0(plotprefix, ".png"),
                width = plotargs$devpars$width,
                height = plotargs$devpars$height,
                res = plotargs$devpars$res
            )
            print(p)
            dev.off()

            reporter$add(
                list(kind = "descr", content = plotargs$descr %||% plot),
                reporter$image(plotprefix, "png", FALSE, kind = "image"),
                h1 = h1,
                h2 = h2
            )
        } else {
            stop("`case` level plot type not supported yet: ", plotargs$plot_type)
        }
    }
}

sapply(names(cases), do_case)

reporter$save(dirname(outdir))
