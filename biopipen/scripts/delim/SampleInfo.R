{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "mutate_helpers.R" | source_r }}

library(rlang)
library(dplyr)

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}
sep <- {{envs.sep | r}}
mutaters <- {{envs.mutaters | r}}
save_mutated <- {{envs.save_mutated | r}}
defaults <- {{envs.defaults | r}}
stats <- {{envs.stats | r}}
exclude_cols <- {{envs.exclude_cols | r}}

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

log_info("Applying mutaters to the data if any ...")
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
add_report(
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
        log_info("- Statistic: {name}")

        case <- cases[[name]]
        info <- casename_info(name, cases, outdir, case_type = "prefix", create = TRUE)
        case <- extract_vars(case, "plot_type", "more_formats", "save_code", "section", subset = "subset_", "devpars", "descr")

        plot_fn <- get_plotthis_fn(plot_type)
        more_formats <- unique(c("png", more_formats))

        if (!is.null(subset_)) {
            case$data <- mutdata %>% dplyr::filter(!!parse_expr(subset_))
        } else {
            case$data <- mutdata
        }

        p <- do_call(plot_fn, case)
        save_plot(p, info$caseprefix, devpars, formats = more_formats)
        if (save_code) {
            save_plotcode(
                p,
                setup = c('library(plotthis)', '', 'load("data.RData")', 'list2env(case, envir = .GlobalEnv)'),
                prefix = info$caseprefix,
                "case"
            )
        }
        report <- list(
            kind = "table_image",
            src = paste0(info$caseprefix, ".png"),
            download = list(),
            name = name,
            descr = descr
        )
        exformats <- setdiff(more_formats, "png")
        if (length(exformats) > 0) {
            report$download <- lapply(exformats, function(fmt) {
                paste0(info$caseprefix, ".", fmt)
            })
        }
        if (save_code) {
            report$download <- c(report$download, list(list(
                src = paste0(info$caseprefix, ".code.zip"),
                tip = "Download the code to reproduce the plot",
                icon = "Code"
            )))
        }

        add_report(report, h1 = "Statistics", ui = "table_of_images:2")
    }
}

save_report(outdir)
