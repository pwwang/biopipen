# Loaded variables: srtfile, outdir, srtobj
library(circlize)

log$info("stats:")

odir <- file.path(outdir, "stats")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)

do_one_stats <- function(name) {
    log$info("- Case: {name}")

    case <- list_update(stats_defaults, stats[[name]])
    extract_vars(case, "devpars", "more_formats", "save_code", "save_data", "subset")

    if (!is.null(subset)) {
        case$object <- srtobj %>% filter(!!parse_expr(subset))
    } else {
        case$object <- srtobj
    }

    info <- case_info(name, odir, is_dir = FALSE, create = TRUE)
    p <- do_call(gglogger::register(CellStatPlot), case)
    save_plot(p, info$prefix, devpars, formats = c("png", more_formats))
    if (save_code) {
        save_plotcode(p, info$prefix,
            setup = c("library(scplotter)", "load('data.RData')", "invisible(list2env('case'))"),
            "case",
            auto_data_setup = FALSE)
    }

    if (save_data) {
        if (!inherits(p$data, "data.frame") && !inherits(p$data, "matrix")) {
            stop("'save_data = TRUE' is not supported for plot_type: ", case$plot_type)
        }
        write.table(p$data, paste0(info$prefix, ".data.txt"), sep="\t", quote=FALSE, row.names=FALSE)
        reporter$add2(
            list(
                name = "Plot",
                contents = list(
                    reporter$image(
                        info$prefix, more_formats, save_code, kind = "image")
                )
            ),
            list(
                name = "Data",
                contents = list(
                    list(
                        kind = "descr",
                        content = "Data used directly for the plot"
                    ),
                    list(
                        kind = "table",
                        src = paste0(info$prefix, ".data.txt"),
                        data = list(nrows = 100)
                    )
                )
            ),
            hs = c(info$section, info$name),
            ui = "tabs"
        )
    } else {
        reporter$add2(
            reporter$image(info$prefix, more_formats, save_code, kind = "image"),
            hs = c(info$section, info$name)
        )
    }
}

sapply(names(stats), do_one_stats)
