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

    return (trimws(unlist(strsplit(features, ","))))
}

do_one_features <- function(name) {
    log$info("- Case: {name}")

    case <- list_update(features_defaults, features[[name]])
    case <- extract_vars(
        case,
        "devpars", "more_formats", "save_code", "save_data", "order_by",
        "subset", "features", "descr")

    if (!is.null(subset)) {
        case$object <- srtobj %>% filter(!!parse_expr(subset))
    } else {
        case$object <- srtobj
    }

    if (exists("order_by") && !is.null(order_by)) {
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
    case$features <- .get_features(features, case$object)
    p <- do_call(gglogger::register(FeatureStatPlot), case)
    save_plot(p, info$prefix, devpars, formats = c("png", more_formats))
    if (save_code) {
        save_plotcode(p, info$prefix,
            setup = c("library(scplotter)", "load('data.RData')", "invisible(list2env('case'))"),
            "case",
            auto_data_setup = FALSE)
    }
    if (exists("descr") && !is.null(descr)) {
        reporter$add2(
            list(
                kind = "descr",
                content = descr
            ),
            hs = c(info$section, info$name)
        )
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

sapply(names(features), do_one_features)
