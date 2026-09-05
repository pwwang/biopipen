library(tidyr)
library(plotthis)
library(biopipen.utils)

covfiles <- {{in.covfiles | each: str | r}}
outdir <- {{out.outdir | str | r}}
joboutdir <- {{job.outdir | r}}

envs <- {{envs | r: todot="-"}}

log <- get_logger()
reporter <- get_reporter()

envs <- extract_vars(envs, "cases", "groups", "save_data", allow_nonexisting = TRUE)

cases <- expand_cases(cases, envs)

log$info("Reading coverage files ...")
covdata <- NULL
for (covfile in covfiles) {
    log$info("- {basename(covfile)} ...")

    tmp <- read.table(covfile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if (ncol(tmp) < 7) {
        stop("Coverage file {covfile} has less than 7 columns. Please check the file format.")
    }
    colnames(tmp)[1:3] <- c("Chrom", "Start", "End")
    colnames(tmp)[ncol(tmp) - 4 + 1:4] <- c("NFeatures", "NRegionBases", "RegionSize", "FracRegionBases")
    tmp$Region <- paste0(tmp$Chrom, ":", tmp$Start, "-", tmp$End)
    tmp$Sample <- sub("\\.coverage$", "", tools::file_path_sans_ext(basename(covfile)))
    if (!is.null(groups) && length(groups) > 0) {
        for (group_name in names(groups)) {
            if (tmp$Sample[1] %in% groups[[group_name]] || paste0(tmp$Sample[1], ".coverage") %in% groups[[group_name]]) {
                tmp$Group <- group_name
                break
            } else {
                tmp$Group <- NA
            }
        }
    } else {
        tmp$Group <- NA
    }
    tmp <- tmp[, c("Region", "Sample", "Group", "NFeatures", "NRegionBases", "RegionSize", "FracRegionBases")]
    covdata <- rbind(covdata, tmp)
}

if (save_data) {
    write.table(covdata, file = file.path(outdir, "coverage.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

    reporter$add(
        list(kind = "table", src = file.path(outdir, "coverage.txt"), data = list(nrows = 100)),
        h1 = "Coverage Data"
    )
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

do_case <- function(name) {
    log$info("Case: {name}")

    info <- case_info(name, outdir, is_dir = FALSE, create = FALSE)
    case <- extract_vars(
        cases[[name]],
        "plot_type", "devpars", "more_formats", "save_code", "descr",
        allow_nonexisting = TRUE
    )

    plot_fn <- get_plotthis_fn(plot_type)

    case$data <- covdata
    p <- do_call(plot_fn,  case)

    save_plot(p, info$prefix, devpars, formats = c("png", more_formats))
    if (save_code) {
        save_plotcode(p, info$prefix,
            setup = c("library(plotthis)", "load('data.RData')", "invisible(list2env(case, envir = .GlobalEnv))"),
            "case",
            auto_data_setup = FALSE)
    }

    reporter$add2(
        list(kind = "descr", content = descr),
        reporter$image(info$prefix, more_formats, save_code, kind = "image"),
        hs = c(info$section, info$name)
    )
}

invisible(sapply(names(cases), do_case))
reporter$save(joboutdir)
