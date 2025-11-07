library(rlang)
library(Seurat)
library(tidyseurat)
library(biopipen.utils)

srtfile <- {{in.srtobj | r}}  # nolint
outdir <- {{out.outdir | r}}  # nolint
joboutdir <- {{job.outdir | r}}  # nolint
mutaters <- {{envs.mutaters | r}}  # nolint
group_by <- {{envs.group_by | default: envs["group-by"] | default: None | r}}  # nolint
ident_1 <- {{envs.ident_1 | default: envs["ident-1"] | default: None | r}}  # nolint
ident_2 <- {{envs.ident_2 | default: envs["ident-2"] | default: None | r}}  # nolint
each <- {{envs.each | r}}  # nolint
subset <- {{envs.subset | r}}  # nolint
gmtfile <- {{envs.gmtfile | r}}  # nolint
method <- {{envs.method | r}}  # nolint
top <- {{envs.top | r}}  # nolint
minsize <- {{envs.minSize | default: envs.minsize | r}}  # nolint
maxsize <- {{envs.maxSize | default: envs.maxsize | r}}  # nolint
eps <- {{envs.eps | r}}  # nolint
alleach_plots_defaults <- {{envs.alleach_plots_defaults | r}}  # nolint
alleach_plots <- {{envs.alleach_plots | r}}  #
ncores <- {{envs.ncores | r}}  # nolint
rest <- {{envs.rest | r: todot="-"}}  # nolint
cases <- {{envs.cases | r: todot="-"}}  # nolint

log <- get_logger()
reporter <- get_reporter()

alleach_plots <- lapply(alleach_plots, function(x) {
    list_update(alleach_plots_defaults, x)
})

log$info("Reading Seurat object ...")
srtobj <- read_obj(srtfile)

if (!is.null(mutaters) && length(mutaters) > 0) {
    log$info("Mutating metadata columns ...")
    srtobj@meta.data <- srtobj@meta.data %>% mutate(!!!lapply(mutaters, parse_expr))
}

defaults <- list(
    group_by = group_by,
    ident_1 = ident_1,
    ident_2 = ident_2,
    each = each,
    subset = subset,
    gmtfile = gmtfile,
    method = method,
    top = top,
    minsize = minsize,
    maxsize = maxsize,
    eps = eps,
    alleach_plots_defaults = alleach_plots_defaults,
    alleach_plots = alleach_plots,
    ncores = ncores,
    rest = rest
)

expand_each <- function(name, case) {
    outcases <- list()

    case$group_by <- case$group_by %||% GetIdentityColumn(srtobj)

    if (is.null(case$each) || is.na(case$each) || nchar(case$each) == 0 || isFALSE(each)) {
        if (length(case$alleach_plots) > 0) {
            stop("Cannot perform `alleach_plots` without `each` defined.")
        }

        outcases[[name]] <- case
    } else {
        eachs <- if (!is.null(case$subset)) {
            srtobj@meta.data %>%
                filter(!!parse_expr(case$subset)) %>%
                pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
        } else {
            srtobj@meta.data %>%
                pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
        }

        if (length(cases) == 0 && name == "GSEA") {
            prefix <- case$each
        } else {
            prefix <- paste0(name, " (", case$each, ")")
        }

        for (each in eachs) {
            newname <- paste0(prefix, "::", each)
            newcase <- case

            newcase$original_case <- paste0(name, " (all ", case$each,")")
            newcase$each_name <- case$each
            newcase$each <- each

            newcase$alleach_plots_defaults <- NULL
            newcase$alleach_plots <- NULL

            if (!is.null(case$subset)) {
                newcase$subset <- paste0(case$subset, " & ", bQuote(case$each), " == '", each, "'")
            } else {
                newcase$subset <- paste0(bQuote(case$each), " == '", each, "'")
            }

            outcases[[newname]] <- newcase
        }

        if (length(case$alleach_plots) > 0) {
            newcase <- case

            newcase$gseas <- list()
            newcase$alleach_plots <- lapply(
                newcase$alleach_plots,
                function(x) { list_update(newcase$alleach_plots_defaults, x) }
            )

            outcases[[paste0(name, " (all ", case$each,")")]] <- newcase
        }
    }
    outcases
}

log$info("Expanding cases...")
cases <- expand_cases(cases, defaults, expand_each, default_case = "GSEA")


ensure_sobj <- function(expr, allow_empty) {
    tryCatch({ expr }, error = function(e) {
        if (allow_empty) {
            log$warn("  Ignoring this case: {e$message}")
            return(NULL)
        } else {
            stop(e)
        }
    })
}

do_case <- function(name) {
    log$info("- Processing case: {name} ...")
    case <- cases[[name]]
    info <- case_info(name, outdir, create = TRUE)

    if (!is.null(case$gseas)) {

        if (length(case$gseas) == 0) {
            log$warn("  No GSEA results found for case {name}. Skipping.")
            return(invisible(NULL))
        }

        each_levels <- names(case$gseas)
        gseas <- do_call(rbind, lapply(each_levels, function(x) {
            gsea_df <- case$gseas[[x]]
            if (nrow(gsea_df) > 0) {
                gsea_df[[case$each]] <- x
            } else {
                gsea_df[[case$each]] <- character(0)  # Empty case
            }
            gsea_df
        }))
        gseas[[case$each]] <- factor(gseas[[case$each]], levels = each_levels)

        for (plotname in names(case$alleach_plots)) {
            plotargs <- case$alleach_plots[[plotname]]
            plotargs <- extract_vars(plotargs, "devpars")
            plotargs$gsea_results <- gseas
            plotargs$group_by <- case$each
            if (plotargs$plot_type == "heatmap") {
                plotargs$show_row_names <- plotargs$show_row_names %||% TRUE
                plotargs$show_column_names <- plotargs$show_column_names %||% TRUE
            }

            p <- do_call(VizGSEA, plotargs)

            outprefix <- file.path(info$prefix, paste0("all.", slugify(plotname)))
            save_plot(p, outprefix, devpars, formats = "png")
            reporter$add2(
                list(kind = "descr", content = paste0("Pathways for all ", case$each, ".")),
                list(kind = "image", src = paste0(outprefix, ".png")),
                hs = c(info$section, info$name),
                hs2 = plotname
            )
        }

        return(invisible(NULL))
    }

    allow_empty = !is.null(case$each)
    # prepare expression matrix
    log$info("  Preparing expression matrix...")
    sobj <- ensure_sobj({ srtobj %>% filter(!is.na(!!sym(case$group_by))) }, allow_empty)
    if (is.null(sobj)) {
        reporter$add2(
            list(
                kind = "error",
                content = paste0("No cells with non-NA `", case$group_by, "` in the Seurat object.")
            ),
            hs = c(info$section, info$name)
        )
        return(NULL)
    }

    if (!is.null(case$subset)) {
        sobj <- ensure_sobj({ sobj %>% filter(!!!parse_exprs(case$subset)) }, allow_empty)
        if (is.null(sobj)) {
            reporter$add2(
                list(
                    kind = "error",
                    content = paste0("No cells with non-NA `", case$group_by, "` in the Seurat object.")
                ),
                hs = c(info$section, info$name)
            )
            return(NULL)
        }
    }
    if (!is.null(case$ident_2)) {
        sobj <- ensure_sobj({ sobj %>% filter(!!sym(case$group_by) %in% c(case$ident_1, case$ident_2)) }, allow_empty)
        if (is.null(sobj)) {
            reporter$add2(
                list(
                    kind = "error",
                    content = paste0("No cells with non-NA `", case$group_by, "` in the Seurat object.")
                ),
                hs = c(info$section, info$name)
            )
            return(NULL)
        }
    }

    allclasses <- sobj@meta.data[, case$group_by, drop = TRUE]
    if (is.null(case$ident_2)) {
        case$ident_2 <- "Other"
        allclasses[allclasses != case$ident_1] <- "Other"
    }
    exprs <- GetAssayData(sobj, layer = "data")

    # get preranks
    log$info("  Getting preranks...")
    ranks <- RunGSEAPreRank(exprs, allclasses, case$ident_1, case$ident_2, case$method)
    write.table(
        as.data.frame(ranks),
        file.path(info$prefix, "fgsea.rank.txt"),
        row.names = TRUE,
        col.names = TRUE,
        sep = "\t",
        quote = FALSE
    )
    if (all(is.na(ranks))) {
        log$warn("  All gene ranks are NA. It's probably due to high missing rate in the data.")
        log$warn("  Case ignored, you may also try a different ranking method.")
        reporter$add2(
            list(
                kind = "error",
                content = "All gene ranks are NA. It's probably due to high missing rate in the data."
            ),
            hs = c(info$section, info$name)
        )
        return(invisible(NULL))
    }

    # run fgsea
    log$info("  Running fgsea...")
    case$rest$ranks <- ranks
    case$rest$genesets <- ParseGMT(case$gmtfile)
    case$rest$minSize <- case$rest$minSize %||% case$rest$minsize %||% case$minsize
    case$rest$maxSize <- case$rest$maxSize %||% case$rest$maxsize %||% case$maxsize
    case$rest$eps <- case$eps
    case$rest$nproc <- case$ncores
    case$rest$minsize <- NULL
    case$rest$maxsize <- NULL
    result <- do_call(RunGSEA, case$rest)
    write.table(
        result,
        file.path(info$prefix, "fgsea.tsv"),
        row.names = FALSE,
        col.names = TRUE,
        sep = "\t",
        quote = FALSE
    )

    aspect.ratio <- sqrt(case$top) / sqrt(10)
    p_summary <- VizGSEA(
        result,
        plot_type = "summary",
        top_term = case$top,
        aspect.ratio = aspect.ratio
    )
    save_plot(
        p_summary,
        file.path(info$prefix, "summary"),
        devpars = list(res = 100, height = attr(p_summary, "height") * 100 / 1.5, width = attr(p_summary, "width") * 100),
        formats = "png"
    )

    p_gsea <- VizGSEA(
        result,
        plot_type = "gsea",
        gs = result$pathway[1:min(case$top, nrow(result))]
    )
    save_plot(
        p_gsea,
        file.path(info$prefix, "pathways"),
        devpars = list(res = 100, height = attr(p_gsea, "height") * 100, width = attr(p_gsea, "width") * 100),
        formats = "png"
    )


    reporter$add2(
        list(
            name = paste0("Table (", case$ident_1, " vs ", case$ident_2, ")"),
            contents = list(
                list(kind = "descr", content = paste0(
                    "Showing top 50 pathways by padj in descending order. ",
                    "Use 'Download the entire data' button to download all pathways."
                )),
                list(kind = "table", src = file.path(info$prefix, "fgsea.tsv"), data = list(nrows = 50))
            )
        ),
        list(
            name = "Summary Plot",
            contents = list(
                list(kind = "descr", content = paste0("Showing top ", case$top, " pathways.")),
                list(kind = "image", src = file.path(info$prefix, "summary.png"))
            )
        ),
        list(
            name = "GSEA Plots",
            contents = list(
                list(kind = "descr", content = paste0("Showing top ", case$top, " pathways.")),
                list(kind = "image", src = file.path(info$prefix, "pathways.png"))
            )
        ),
        hs = c(info$section, info$name),
        ui = "tabs"
    )

    if (!is.null(case$original_case) && !is.null(cases[[case$original_case]])) {
        cases[[case$original_case]]$gseas[[case$each]] <<- result
    }

    invisible()
}

sapply(names(cases), function(name) do_case(name))

reporter$save(joboutdir)
