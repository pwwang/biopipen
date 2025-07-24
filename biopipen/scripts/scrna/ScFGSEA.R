library(rlang)
library(Seurat)
library(tidyseurat)
library(biopipen.utils)

srtfile <- {{in.srtobj | r}}  # nolint
outdir <- {{out.outdir | r}}  # nolint
joboutdir <- {{job.outdir | r}}  # nolint
mutaters <- {{envs.mutaters | r}}  # nolint
group.by <- {{envs["group-by"] | r}}  # nolint
ident.1 <- {{envs["ident-1"] | r}}  # nolint
ident.2 <- {{envs["ident-2"] | r}}  # nolint
each <- {{envs.each | r}}  # nolint
subset <- {{envs.subset | r}}  # nolint
gmtfile <- {{envs.gmtfile | r}}  # nolint
method <- {{envs.method | r}}  # nolint
top <- {{envs.top | r}}  # nolint
minsize <- {{envs.minSize | default: envs.minsize | r}}  # nolint
maxsize <- {{envs.maxSize | default: envs.maxsize | r}}  # nolint
eps <- {{envs.eps | r}}  # nolint
ncores <- {{envs.ncores | r}}  # nolint
rest <- {{envs.rest | r: todot="-"}}  # nolint
cases <- {{envs.cases | r: todot="-"}}  # nolint

log <- get_logger()
reporter <- get_reporter()

log$info("Reading Seurat object ...")
srtobj <- read_obj(srtfile)
if (!"Identity" %in% colnames(srtobj@meta.data)) {
    srtobj@meta.data$Identity <- Idents(srtobj)
}

if (!is.null(mutaters) && length(mutaters) > 0) {
    log$info("Mutating metadata columns ...")
    srtobj@meta.data <- srtobj@meta.data %>% mutate(!!!lapply(mutaters, parse_expr))
}

defaults <- list(
    group.by = group.by,
    ident.1 = ident.1,
    ident.2 = ident.2,
    each = each,
    subset = subset,
    gmtfile = gmtfile,
    method = method,
    top = top,
    minsize = minsize,
    maxsize = maxsize,
    eps = eps,
    ncores = ncores,
    rest = rest
)

expand_each <- function(name, case) {
    outcases <- list()

    case$group.by <- case$group.by %||% "Identity"

    if (is.null(case$each) || is.na(case$each) || nchar(case$each) == 0 || isFALSE(each)) {
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
            name <- case$each
        }

        for (each in eachs) {
            newname <- paste0(case$each, "::", each)
            newcase <- case

            newcase$original_case <- name
            newcase$each_name <- case$each
            newcase$each <- each

            if (!is.null(case$subset)) {
                newcase$subset <- paste0(case$subset, " & ", bQuote(case$each), " == '", each, "'")
            } else {
                newcase$subset <- paste0(bQuote(case$each), " == '", each, "'")
            }

            outcases[[newname]] <- newcase
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

    allow_empty = !is.null(case$each)
    # prepare expression matrix
    log$info("  Preparing expression matrix...")
    sobj <- ensure_sobj({ srtobj %>% filter(!is.na(!!sym(case$group.by))) }, allow_empty)
    if (is.null(sobj)) {
        reporter$add2(
            list(
                kind = "error",
                content = paste0("No cells with non-NA `", case$group.by, "` in the Seurat object.")
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
                    content = paste0("No cells with non-NA `", case$group.by, "` in the Seurat object.")
                ),
                hs = c(info$section, info$name)
            )
            return(NULL)
        }
    }
    if (!is.null(case$ident.2)) {
        sobj <- ensure_sobj({ sobj %>% filter(!!sym(case$group.by) %in% c(case$ident.1, case$ident.2)) }, allow_empty)
        if (is.null(sobj)) {
            reporter$add2(
                list(
                    kind = "error",
                    content = paste0("No cells with non-NA `", case$group.by, "` in the Seurat object.")
                ),
                hs = c(info$section, info$name)
            )
            return(NULL)
        }
    }

    allclasses <- sobj@meta.data[, case$group.by, drop = TRUE]
    if (is.null(case$ident.2)) {
        case$ident.2 <- "Other"
        allclasses[allclasses != case$ident.1] <- "Other"
    }
    exprs <- GetAssayData(sobj, layer = "data")

    # get preranks
    log$info("  Getting preranks...")
    ranks <- RunGSEAPreRank(exprs, allclasses, case$ident.1, case$ident.2, case$method)
    write.table(
        ranks,
        file.path(info$prefix, "fgsea.rank"),
        row.names = FALSE,
        col.names = TRUE,
        sep = "\t",
        quote = FALSE
    )
    if (all(is.na(ranks))) {
        if (length(allclasses) < 100) {
            log$warn("  Ignoring this case because all gene ranks are NA and there are <100 cells.")
            reporter$add2(
                list(
                    kind = "error",
                    content = paste0("Not enough cells (n = ", length(allclasses), ") to run fgsea.")
                ),
                hs = c(info$section, info$name)
            )
            return(NULL)
        } else {
            stop(paste0(
                "All gene ranks are NA (# cells = ",
                length(allclasses),
                "). ",
                "It's probably due to high missing rate in the data. ",
                "You may want to try a different `envs$method` for pre-ranking."
            ))
        }
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

    p_summary <- VizGSEA(
        result,
        plot_type = "summary",
        top_term = case$top
    )
    save_plot(
        p_summary,
        file.path(info$prefix, "summary"),
        devpars = list(res = 100, height = attr(p_summary, "height") * 100, width = attr(p_summary, "width") * 100),
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
            name = "Table",
            contents = list(
                list(kind = "descr", content = paste0(
                    "Showing top 50 pathways by padj in descending order. ",
                    "Use 'Download the entire data' button to download all pathways."
                )),
                list(kind = "table", src = file.path(info$prefix, "fgsea"), data = list(nrows = 50))
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
}

sapply(sort(names(cases)), function(name) do_case(name))

reporter$save(joboutdir)
