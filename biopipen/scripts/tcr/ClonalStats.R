{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "scripts", "tcr", "ScRep-common.R" | source_r }}

library(scplotter)

screpfile <- {{in.screpfile | quote}}
outdir <- {{out.outdir | quote}}
joboutdir <- {{job.outdir | quote}}
envs <- {{envs | r: todot="-"}}
mutaters <- envs$mutaters
cases <- envs$cases
envs$mutaters <- NULL
envs$cases <- NULL

VIZ_TYPES <- list(
    volume = "Number of Clones",
    abundance = "Clonal Abundance",
    length = "Clonal Sequence Length",
    residency = "Clonal Residency",
    composition = "Clonal Composition",
    overlap = "Clonal Overlap",
    diversity = "Clonal Diversity",
    geneusage = "Gene Usage",
    positional = "Positional Properties",
    kmer = "Kmer Analysis",
    rarefaction = "Rarefaction Analysis"
)

log_info("Loading scRepertoire object ...")
screp <- readRDS(screpfile)

log_info("Applying mutaters if any ...")
screp <- screp_mutate(screp, mutaters)

log_info("Making cases ...")
cases <- expand_cases(cases, envs)
viz_types <- list()
for (name in names(cases)) {
    case <- cases[[name]]
    if (is.null(case$viz_type)) {
        stop("Error: Visualization type is not defined for case '", name, "'")
    }
    if (!case$viz_type %in% names(VIZ_TYPES)) {
        stop("Error: Unknown visualization type '", case$viz_type, "' for case '", name,
             "'. Available types: ", paste(names(VIZ_TYPES), collapse = ", "))
    }
    if (is.null(case$section)) {
        viz_types[[case$viz_type]] <- viz_types[[case$viz_type]] %||% 0
        viz_types[[case$viz_type]] <- viz_types[[case$viz_type]] + 1
    }
}
for (name in names(cases)) {
    case <- cases[[name]]
    if (is.null(case$section) && viz_types[[case$viz_type]] > 1) {
        cases[[name]]$section <- VIZ_TYPES[[case$viz_type]]
    }
}

do_case <- function(name, case) {
    log_info("- Processing case: {name}")
    section <- case$section
    info <- casename_info(name, cases, outdir, section = section, case_type = "prefix", create = TRUE)

    case <- extract_vars(case, "viz_type", "devpars", "more_formats", "save_code", "section", subset = "subset_")

    if (!is.null(subset_)) {
        case$data <- screp_subset(screp, subset_)
    } else {
        case$data <- screp
    }

    plot_fn <- paste0("Clonal", tools::toTitleCase(viz_type), "Plot")
    plot_fn <- utils::getFromNamespace(plot_fn, "scplotter")
    if (is.null(plot_fn)) {
        stop("Error: Unknown visualization type: ", viz_type)
    }

    p <- do_call(plot_fn, case)
    save_plot(p, info$caseprefix, devpars, formats = unique(c("png", more_formats)))

    report <- list(
        kind = "table_image",
        src = paste0(info$caseprefix, ".png"),
        download = list(),
        name = name
    )
    exformats <- setdiff(more_formats, "png")
    if (length(exformats) > 0) {
        report$download <- lapply(exformats, function(fmt) {
            paste0(info$caseprefix, ".", fmt)
        })
    }

    if (isTRUE(save_code)) {
        save_plotcode(
            p,
            setup = c('library(scplotter)', '', 'load("data.RData")'),
            prefix = info$caseprefix,
            "case"
        )
        report$download <- c(report$download, list(list(
            src = paste0(info$caseprefix, ".code.zip"),
            tip = "Download the code to reproduce the plot",
            icon = "Code"
        )))
    }

    if (is.null(section)) {
        h1 <- html_escape(name)
        h2 <- "#"
    } else {
        h1 <- info$h1
        h2 <- info$h2
    }
    add_report(report, h1 = h1, h2 = h2, ui = "table_of_images:2")
}

lapply(names(cases), function(name) do_case(name, cases[[name]]))

save_report(joboutdir)
