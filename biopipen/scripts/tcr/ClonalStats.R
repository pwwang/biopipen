{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "scripts", "tcr", "ScRep-common.R" | source_r }}

library(scplotter)

screpfile <- {{in.screpfile | quote}}
outdir <- {{out.outdir | quote}}
envs <- {{envs | r}}
mutaters <- envs$mutaters
cases <- envs$cases
envs$mutaters <- NULL
envs$cases <- NULL

log_info("Loading scRepertoire object ...")
screp <- readRDS(screpfile)

log_info("Applying mutaters if any ...")
screp <- screp_mutate(screp, mutaters)

log_info("Making cases ...")
cases <- expand_cases(cases, envs)

do_case <- function(name, case) {
    log_info("- Processing case: {name}")
    info <- casename_info(name, cases, outdir, case_type = "prefix", create = TRUE)

    viz_type <- case$viz_type
    if (is.null(viz_type)) {
        stop("Error: Visualization type is not defined for case: ", case$name)
    }

    devpars <- case$devpars %||% list()
    case$devpars <- NULL
    subset_ <- case$subset
    case$subset <- NULL
    pdf_ <- case$pdf
    case$pdf <- NULL
    code <- case$code
    case$code <- NULL
    case$section <- NULL

    if (!is.null(subset_)) {
        case$data <- screp_subset(screp, subset_)
    } else {
        case$data <- screp
    }

    plot_fn <- paste0("Clonal", tools::toTitleCase(case$viz_type), "Plot")
    plot_fn <- utils::getFromNamespace(plot_fn, "scplotter")
    if (is.null(plot_fn)) {
        stop("Error: Unknown visualization type: ", case$viz_type)
    }

    p <- do_call(plot_fn, case)
    save_plot(p, info$caseprefix, devpars, formats = if (pdf_) c("png", "pdf") else "png")

    if (isTRUE(code)) {
        save_plotcode(
            p,
            setup = c('library(scplotter)', '', 'load("data.RData")'),
            prefix = info$caseprefix,
            "case"
        )
    }
}

lapply(names(cases), function(name) do_case(name, cases[[name]]))
