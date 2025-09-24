library(rlang)
library(dplyr)
library(scplotter)
library(biopipen.utils)

cccfile <- {{ in.cccfile | r }}
outdir <- {{ out.outdir | r }}
joboutdir <- {{ job.outdir | r }}
envs <- {{ envs | r }}
envs <- extract_vars(
    envs,
    "magnitude", "specificity", "devpars", "subset", "cases", "more_formats", "descr"
)

ccc <- read.table(cccfile, header=TRUE, sep="\t", check.names = FALSE)

if (length(ccc) == 0) {
    stop("No data found in the input file: ", cccfile)
}

defaults <- list(
    magnitude = NULL,
    specificity = NULL,
    subset = subset,
    descr = descr,
    more_formats = more_formats,
    devpars = list(res = 100)
)

cases <- expand_cases(cases, defaults, default_case = "Cell-Cell Communication")
log <- get_logger()
reporter <- get_reporter()

do_case <- function(name) {
    log$info("- Case: {name}")
    case <- cases[[name]]
    info <- case_info(name, outdir, is_dir = FALSE)
    case <- extract_vars(case, subset_ = "subset", "devpars", "more_formats", "descr")

    case$data <- ccc
    if (!is.null(subset_)) {
        case$data <- ccc %>% dplyr::filter(!!parse_expr(subset_))
    }

    if (identical(case$plot_type, "table")) {
        write.table(
            case$data,
            file = paste0(info$prefix, ".txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )
        report <- list(
            kind = "table",
            data = list(nrows = 100),
            src = paste0(info$prefix, ".txt")
        )
        reporter$add2(report, hs = c(info$section, info$name))
        return()
    }

    if (is.null(case$magnitude)) {
        case$magnitude <- NULL
    }
    if (is.null(case$specificity)) {
        case$specificity <- NULL
    }
    p <- do_call(scplotter::CCCPlot, case)
    save_plot(
        p, info$prefix,
        devpars = devpars, formats = unique(c("png", more_formats))
    )

    report <- list(
        kind = "table_image",
        src = paste0(info$prefix, ".png"),
        download = list(),
        descr = html_escape(descr),
        name = html_escape(info$name)
    )
    exformats <- setdiff(more_formats, "png")
    if (length(exformats) > 0) {
        report$download <- lapply(exformats, function(fmt) {
            paste0(info$prefix, ".", fmt)
        })
    }
    reporter$add2(report, hs = c(info$section, info$name), ui = "table_of_images:2")
}

sapply(names(cases), do_case)

reporter$save(joboutdir)
