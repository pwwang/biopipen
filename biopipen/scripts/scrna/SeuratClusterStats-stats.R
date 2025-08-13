# Loaded variables: srtfile, outdir, srtobj

log$info("stats:")

odir <- file.path(outdir, "stats")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)



do_one_stats <- function(name) {
    log$info("- Case: {name}")

    case <- list_update(stats_defaults, stats[[name]])
    case <- extract_vars(case, "devpars", "more_formats", "save_code", "save_data", "subset", "descr")

    if (!is.null(subset)) {
        case$object <- srtobj %>% filter(!!parse_expr(subset))
    } else {
        case$object <- srtobj
    }
    ident <- case$ident %||% GetIdentityColumn(case$object)
    groupings <- unique(c(case$group_by, case$rows_by, case$columns_by, case$pie_group_by, ident))
    if (length(groupings) > 0) {
        for (g in groupings) {
            case$object <- filter(case$object, !is.na(!!sym(g)))
        }
    }

    info <- case_info(name, odir, is_dir = FALSE, create = TRUE)
    p <- do_call(gglogger::register(CellStatPlot), case)
    save_plot(p, info$prefix, devpars, formats = c("png", more_formats))
    if (save_code) {
        save_plotcode(p, info$prefix,
            setup = c("library(scplotter)", "load('data.RData')", "invisible(list2env(case, envir = .GlobalEnv))"),
            "case",
            auto_data_setup = FALSE)
    }

    frac <- case$frac %||% "none"
    default_descr <- glue(
        "The {case$plot_type} plot shows the distribution of cells across categories defined by '{ident}'",
        "{if (!is.null(case$group_by)) glue(', grouped by {case$group_by}') else ''}",
        "{if (!is.null(case$split_by)) glue(', and split by {case$split_by}') else ''}. ",
        "The values represent ",
        "{if (frac == 'none') 'the number of cells' else glue('the fraction of cells calculated by \"{frac}\"')}. "
    )
    if (!is.null(case$comparisons)) {
        default_descr <- paste0(
            default_descr,
            glue("Statistical comparisons were performed between groups using \"{case$pairwise_method %||% 'wilcox.test'}\" method.")
        )
    }
    if (save_data) {
        pdata <- attr(p, "data") %||% p$data
        if (!inherits(pdata, "data.frame") && !inherits(pdata, "matrix")) {
            stop("'save_data = TRUE' is not supported for plot_type: ", case$plot_type)
        }
        write.table(pdata, paste0(info$prefix, ".data.txt"), sep="\t", quote=FALSE, row.names=FALSE)
        reporter$add2(
            list(
                name = "Plot",
                contents = list(
                    list(
                        kind = "descr",
                        content = case$descr %||% default_descr
                    ),
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
            list(kind = "descr", content = case$descr %||% default_descr),
            reporter$image(info$prefix, more_formats, save_code, kind = "image"),
            hs = c(info$section, info$name)
        )
    }
}

sapply(names(stats), do_one_stats)
