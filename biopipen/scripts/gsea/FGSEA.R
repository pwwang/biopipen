library(rlang)
library(biopipen.utils)

# input & output
infile = {{in.infile | r}}
metafile = {{in.metafile | r}}
outdir = {{out.outdir | r}}
joboutdir = {{job.outdir | r}}

# envs
ncores = {{envs.ncores | r}}
case = {{envs.case | r}}
control = {{envs.control | r}}
gmtfile = {{envs.gmtfile | r}}
method = {{envs.method | r}}
clscol = {{envs.clscol | r}}
top = {{envs.top | r}}
eps = {{envs.eps | r}}
minsize = {{envs.minSize | default: envs.minsize | r}}
maxsize = {{envs.maxSize | default: envs.maxsize | r}}
rest = {{envs.rest | r}}
cases = {{envs.cases | r}}

log <- get_logger()
reporter <- get_reporter()

defaults <- list(
    case = case,
    control = control,
    gmtfile = gmtfile,
    method = method,
    clscol = clscol,
    top = top,
    eps = eps,
    minsize = minsize,
    maxsize = maxsize,
    rest = rest
)
cases <- expand_cases(cases, defaults, default_case = "GSEA")

log$info("Reading input file ...")
indata <- read.table(infile, header=TRUE, stringsAsFactors=FALSE, row.names=1, sep="\t", quote="", check.names=FALSE)

if (!is.null(metafile)) {
    log$info("Reading metadata file ...")
    metadata <- read.table(metafile, header=TRUE, stringsAsFactors=FALSE, row.names=NULL, sep="\t", quote="", check.names=FALSE)
} else {
    metadata <- NULL
}

do_case <- function(name) {
    log$info("Processing case: {name} ...")
    case <- cases[[name]]
    info <- case_info(name, outdir, create = TRUE)

    if (is.null(case$case) && is.null(case$control)) {
        stop("Either `case` or `control` must be specified in the case.")
    }
    if (is.null(case$gmtfile)) {
        stop("`gmtfile` must be specified in the case.")
    }
    if (is.null(case$clscol)) {
        stop("`clscol` must be specified in the case.")
    }
    if (!is.null(metadata) && length(case$clscol) > 1) {
        stop("When `in.metafile` is specified, `envs.clscol` must be a single column name.")
    }
    if (!is.null(metadata)) {
        samples <- colnames(indata)
        if (!"Sample" %in% colnames(metadata)) {
            colnames(metadata)[1] <- "Sample"
        }
        metadata <- metadata[match(samples, metadata$Sample), , drop=FALSE]
        case$clscol <- as.character(metadata[[case$clscol]])
    }
    if (length(unique(case$clscol)) < 2) {
        stop("The `clscol` must have at least two unique values.")
    }
    if (length(unique(case$clscol)) == 2) {
        case$case <- case$case %||% setdiff(unique(case$clscol), case$control)
        case$control <- case$control %||% setdiff(unique(case$clscol), case$case)
    } else {
        if (is.null(case$case) || is.null(case$control)) {
            stop("When `clscol` has more than two unique values, both `case` and `control` must be specified.")
        }
    }
    log$info("- Running pre-ranking ...")
    ranks <- RunGSEAPreRank(
        indata,
        classes = case$clscol,
        case = case$case,
        control = case$control,
        method = case$method
    )
    if (all(is.na(ranks))) {
        if (length(case$clscol) < 10) {
            log$warn("  Ignoring this case because all gene ranks are NA and there are <10 samples.")
            reporter$add2(
                list(
                    kind = "error",
                    content = paste0("Not enough samples (n = ", length(case$clscol), ") to run fgsea.")
                ),
                hs = c(info$section, info$name)
            )
            return(NULL)
        } else {
            stop(paste0(
                "All gene ranks are NA (# samples = ",
                length(case$clscol),
                "). ",
                "It's probably due to high missing rate in the data. ",
                "You may want to try a different `envs$method` for pre-ranking."
            ))
        }
    }

    log$info("- Running GSEA ...")
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

sapply(names(cases), do_case)
reporter$save(joboutdir)
