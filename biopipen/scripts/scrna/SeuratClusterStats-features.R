# Loaded variables: srtfile, outdir, srtobj

# features_defaults = {{envs.features_defaults | r: todot="-"}}
# features = {{envs.features | r: todot="-", skip=1}}
log$info("features:")

odir = file.path(outdir, "features")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)

# highly variable features
hvf <- NULL

.get_features = function(features, object) {
    if (is.null(features)) { features = 20 }
    if (is.numeric(features)) {
        if (!is.null(hvf)) {
            return(hvf[1:features])
        }
        vf <- VariableFeatures(object)
        if (length(vf) == 0) {
            if (DefaultAssay(object) == "SCT") {
                # Still use RNA assay to find variable features
                # See
                # https://github.com/satijalab/seurat/issues/6064
                # https://github.com/satijalab/seurat/issues/8238
                # https://github.com/satijalab/seurat/issues/5761
                vf <- FindVariableFeatures(object, nfeatures = features, assay = "RNA")
            } else {
                vf <- FindVariableFeatures(object, nfeatures = features)
            }
        }
        hvf <<- vf
        return(hvf[1:features])
    }
    if (is.character(features) && length(features) > 1) {
        return (features)
    }
    if (is.character(features) && startsWith(features, "file://")) {
        return (read.table(
            substring(features, 8),
            sep = "\t",
            header = FALSE,
            row.names = NULL,
            check.names = FALSE
        )$V1)
    }

    if (is.null(features)) {
        if (is.null(default_features)) {
            return (default[1:20])
        } else {
            return (default_features)
        }
    }

    if (is.list(features)) {
        return(lapply(features, function(x) {.get_features(x, object) }))
    }

    return (trimws(unlist(strsplit(features, ","))))
}

do_one_features <- function(name) {
    log$info("- Case: {name}")

    case <- list_update(features_defaults, features[[name]])
    case <- extract_vars(
        case,
        "devpars", "more_formats", "save_code", "save_data", "order_by",
        "subset", "features", "descr",
        allow_nonexisting = TRUE)

    if (!is.null(subset)) {
        case$object <- srtobj %>% filter(!!parse_expr(subset))
    } else {
        case$object <- srtobj
    }

    if (exists("order_by") && !is.null(order_by)) {
        case$ident <- case$ident %||% GetIdentityColumn(case$object)
        if (length(order_by) < 2) {
            clusters <- case$object@meta.data %>%
                group_by(!!sym(case$ident)) %>%
                arrange(!!parse_expr(order_by)) %>%
                ungroup() %>%
                pull(!!sym(case$ident)) %>%
                unique()

            case$object@meta.data[[case$ident]] <- factor(case$object@meta.data[[case$ident]], levels = clusters)
        } else {
            case$object@meta.data[[case$ident]] <- fct_relevel(case$object@meta.data[[case$ident]], order_by)
        }
    }

    info <- case_info(name, odir, is_dir = FALSE, create = TRUE)

    caching <- Cache$new(
        c(case, list(devpars, more_formats, save_code, save_data, order_by, subset, features, descr)),
        prefix = "biopipen.scrna.SeuratClusterStats.features",
        cache_dir = cache,
        kind = "prefix",
        path = info$prefix
    )

    if (caching$is_cached()) {
        log$info("  plots are cached, restoring ...")
        caching$restore()
    } else {
        case$features <- .get_features(features, case$object)
        p <- tryCatch({
            do_call(gglogger::register(FeatureStatPlot), case)
        }, error = function(e) {
            if (save_code) { stop(e) }
            do_call(FeatureStatPlot, case)
        })
        save_plot(p, info$prefix, devpars, formats = c("png", more_formats))
        if (save_code) {
            save_plotcode(p, info$prefix,
                setup = c("library(scplotter)", "load('data.RData')", "invisible(list2env(case, envir = .GlobalEnv))"),
                "case",
                auto_data_setup = FALSE)
        }

        if (save_data) {
            pdata <- attr(p, "data") %||% p$data
            if (!inherits(pdata, "data.frame") && !inherits(pdata, "matrix")) {
                stop("'save_data = TRUE' is not supported for plot_type: ", case$plot_type)
            }
            write.table(pdata, paste0(info$prefix, ".data.txt"), sep="\t", quote=FALSE, row.names=FALSE)
        }

        caching$save(info$prefix)
    }
    # add reports
    default_descr <- glue(
        "The plot shows the distribution or pattern of the specified features ({paste(case$features %||% features, collapse = ', ')}) ",
        "across cells",
        "{if (!is.null(case$ident)) glue(', identified by \"{case$ident}\"') else ''}",
        "{if (!is.null(case$group_by)) glue(', grouped by \"{case$group_by}\"') else ''}",
        "{if (!is.null(case$split_by)) glue(', and split by \"{case$split_by}\"') else ''}. ",
        "The plot type is '{case$plot_type}', ",
        "{if (case$plot_type == 'dim') 'displaying the features on a dimensional reduction embedding' ",
        " else if (case$plot_type == 'heatmap') 'arranged as a heatmap by rows_name and other grouping variables' ",
        " else if (case$plot_type %in% c('violin', 'box', 'ridge')) 'showing the distribution of feature values by the grouping variables' ",
        " else if (case$plot_type == 'cor') 'showing the correlation between features' ",
        " else 'showing aggregated feature values by the grouping variables'}. ",
        "{if (!is.null(case$facet_by)) glue('Plots are further faceted by \"{case$facet_by}\". ') else ''}",
        "{if (case$plot_type == 'dim') glue('The reduction used is \"{if (!is.null(case$reduction)) case$reduction else DefaultDimReduc(case$object)}\"') else ''}",
        "{if (case$plot_type == 'dim' && !is.null(case$graph)) glue(', with graph \"{case$graph}\" drawn to show cell neighbor edges') else ''}",
        "{if (case$plot_type == 'dim' && !is.null(case$bg_cutoff) && case$bg_cutoff > 0) glue(', and a background cutoff of {case$bg_cutoff}') else ''}",
        "{if (case$plot_type == 'dim') glue(', using dimensions {paste(case$dims %||% 1:2, collapse = \",\")}') else ''}"
    )
    if (!is.null(case$comparisons)) {
        default_descr <- paste0(
            default_descr,
            glue("Statistical comparisons were performed between groups using \"{case$pairwise_method %||% 'wilcox.test'}\" method.")
        )
    }
    reporter$add2(
        list(kind = "descr", content = descr %||% default_descr),
        hs = c(info$section, info$name)
    )

    if (save_data) {
        reporter$add2(
            list(
                name = "Plot",
                contents = list(reporter$image(info$prefix, more_formats, save_code, kind = "image"))
            ),
            list(
                name = "Data",
                contents = list(
                    list(kind = "descr", content = "Data used directly for the plot"),
                    list(kind = "table", src = paste0(info$prefix, ".data.txt"), data = list(nrows = 100))
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

sapply(names(features), do_one_features)
