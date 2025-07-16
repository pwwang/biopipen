# Loaded variables: srtfile, outdir, srtobj

# ngenes_defaults <- {{envs.ngenes_defaults | r: todot="-"}}
# ngenes <- {{envs.ngenes | r: todot="-", skip=1}}
log$info("ngenes:")

odir <- file.path(outdir, "ngenes")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)

do_one_ngenes <- function(name) {
    log$info("- Case: {name}")

    case <- list_update(ngenes_defaults, ngenes[[name]])
    case$devpars <- list_update(ngenes_defaults$devpars, case$devpars)
    case$more_formats <- case$more_formats %||% character(0)
    case$save_code <- case$save_code %||% FALSE
    case$descr <- case$descr %||% name
    case$save_data <- case$save_data %||% FALSE
    case$ylab <- case$ylab %||% "Number of expressed genes"
    case$features <- "Number of expressed genes"
    extract_vars(case, "devpars", "more_formats", "descr", "save_code", "save_data", subset_ = "subset")

    if (!is.null(case$subset)) {
        case$object <- srtobj %>% filter(!!rlang::parse_expr(subset_))
    } else {
        case$object <- srtobj
    }
    case$object <- AddMetaData(case$object, Matrix::colSums(GetAssayData(case$object) > 0), col.name = "Number of expressed genes")

    info <- case_info(name, odir, is_dir = FALSE, create = TRUE)
    p <- do_call(gglogger::register(FeatureStatPlot), case)
    save_plot(p, info$prefix, case$devpars, formats = c("png", more_formats))
    if (save_code) {
        save_plotcode(p, info$prefix,
            setup = c("library(scplotter)", "load('data.RData')", "invisible(list2env(case, envir = .GlobalEnv))"),
            "case",
            auto_data_setup = FALSE
        )
    }
    if (save_data) {
        pdata <- attr(p, "data") %||% p$data
        if (!inherits(pdata, "data.frame") && !inherits(pdata, "matrix")) {
            stop("'save_data = TRUE' is not supported for plot_type: ", case$plot_type)
        }
        write.table(pdata, paste0(info$prefix, ".data.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
        reporter$add2(
            list(
                name = "Plot",
                contents = list(
                    list(kind = "descr", content = case$descr),
                    reporter$image(info$prefix, more_formats, save_code, kind = "image")
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
    }
    else {
        reporter$add2(
            list(kind = "descr", content = case$descr),
            reporter$image(info$prefix, more_formats, save_code, kind = "image"),
            hs = c(info$section, info$name)
        )
    }
}

sapply(names(ngenes), do_one_ngenes)
