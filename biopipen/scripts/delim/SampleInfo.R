library(rlang)
library(dplyr)
library(gglogger)
library(biopipen.utils)
library(plotthis)

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}
joboutdir <- {{job.outdir | r}}
sep <- {{envs.sep | r}}
mutaters <- {{envs.mutaters | r}}
save_mutated <- {{envs.save_mutated | r}}
defaults <- {{envs.defaults | r}}
stats <- {{envs.stats | r}}
exclude_cols <- {{envs.exclude_cols | r}}

log <- get_logger()
reporter <- get_reporter()

if (is.null(exclude_cols)) {
    exclude_cols <- c()
} else {
    exclude_cols <- trimws(unlist(strsplit(exclude_cols, ",")))
}

outdir <- dirname(outfile)
indata <- read.delim(infile, sep = sep, header = TRUE, row.names = NULL)

if (colnames(indata)[1] == "row.names") {
    stop("Wrong number of column names. Do you have the right `sep`?")
}

#' Get plotthis function from plot_type
#'
#' @param plot_type The plot type
#' @param gglogger_register Register the plotthis function to gglogger
#' @param return_name Return the name of the function instead of the function
#' @return The plotthis function
#' @export
get_plotthis_fn <- function(plot_type, gglogger_register = TRUE, return_name = FALSE) {
    fn_name <- switch(plot_type,
        hist = "Histogram",
        histo = "Histogram",
        histogram = "Histogram",
        featuredim = "FeatureDimPlot",
        splitbar = "SplitBarPlot",
        enrichmap = "EnrichMap",
        enrichnet = "EnrichNetwork",
        enrichnetwork = "EnrichNetwork",
        gsea = "GSEAPlot",
        gseasummary = "GSEASummaryPlot",
        gseasum = "GSEASummaryPlot",
        heatmap = "Heatmap",
        network = "Network",
        pie = "PieChart",
        wordcloud = "WordCloudPlot",
        venn = "VennDiagram",
        {
            title_case_plot_type <- tools::toTitleCase(plot_type)
            if (endsWith(title_case_plot_type, "Plot")) {
                title_case_plot_type
            } else if (endsWith(title_case_plot_type, "plot")) {
                paste0(substr(title_case_plot_type, 1, nchar(title_case_plot_type) - 4), "Plot")
            } else {
                paste0(title_case_plot_type, "Plot")
            }
        }
    )
    if (return_name) {
        return(fn_name)
    }
    fn <- tryCatch({
        utils::getFromNamespace(fn_name, "plotthis")
    }, error = function(e) {
        stop("Unknown plot type: ", plot_type)
    })

    if (gglogger_register) {
        gglogger::register(fn, fn_name)
    } else {
        fn
    }
}

log$info("Applying mutaters to the data if any ...")
if (!is.null(mutaters) && length(mutaters) > 0) {
    mutdata <- indata %>%
        mutate(!!!lapply(mutaters, parse_expr))
} else {
    mutdata <- indata
}

write.table(
    if (save_mutated) mutdata else indata,
    file = outfile,
    sep = sep,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)


reporter$add(
    list(
        kind = "descr",
        content = "The samples used in the analysis. Each row is a sample, and columns are the meta information about the sample. This is literally the input sample information file, but the paths to the scRNA-seq and scTCR-seq data are hidden.",
        once = TRUE
    ),
    list(
        kind = "table",
        pageSize = 50,
        data = list(file = outfile, sep = sep, excluded = exclude_cols),
        src = FALSE
    ),
    h1 = "Sample Information"
)

if (length(stats) > 0) {
    cases <- expand_cases(stats, defaults)
    for (name in names(cases)) {
        log$info("- Statistic: {name}")

        case <- cases[[name]]
        info <- case_info(name, outdir, is_dir = FALSE, create = TRUE)
        case <- extract_vars(case, "plot_type", "more_formats", "save_code", "section", "subset", "devpars", "descr")

        plot_fn <- get_plotthis_fn(plot_type)
        more_formats <- unique(c("png", more_formats))

        if (!is.null(subset)) {
            case$data <- mutdata %>% dplyr::filter(!!parse_expr(subset))
        } else {
            case$data <- mutdata
        }

        p <- do_call(plot_fn, case)
        save_plot(p, info$prefix, devpars, formats = more_formats)
        if (save_code) {
            save_plotcode(
                p,
                setup = c('library(plotthis)', '', 'load("data.RData")', 'list2env(case, envir = .GlobalEnv)'),
                prefix = info$caseprefix,
                "case",
                auto_data_setup = FALSE
            )
        }

        reporter$add(
            reporter$image(
                info$prefix,
                c("png", more_formats),
                save_code,
                kind = "table_image"
            ),
            h1 = "Statistics", ui = "table_of_images:2"
        )
    }
}

reporter$save(joboutdir)
