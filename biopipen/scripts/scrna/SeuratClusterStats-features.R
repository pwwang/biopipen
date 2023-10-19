# Loaded variables: srtfile, outdir, srtobj

features_defaults = {{envs.features_defaults | r: todot="-"}}
features = {{envs.features | r: todot="-", skip=1}}

odir = file.path(outdir, "features")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)
report_toc_file = file.path(odir, "report_toc.json")
# Section => list(
#   list(name?, kind, file),
#   ...
# )
report_toc = list()

.add_toc = function(section, toc) {
    if (section %in% names(report_toc)) {
        report_toc[[section]][length(report_toc[[section]])] <<- toc
    } else {
        report_toc[[section]] <<- list(toc)
    }
}

.get_features = function(features) {
    if (is.null(features)) { features = 20 }
    if (is.numeric(features)) {
        return (VariableFeatures(srtobj)[1:features])
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

do_one_features = function(name) {
    print(paste0("Doing features for: ", name))

    case = list_update(features_defaults, features[[name]])
    toc = list()
    if (!is.null(case$section)) { toc$name = name }
    case$devpars = list_update(features_defaults$devpars, features[[name]]$devpars)
    excluded_args = c(
        "section",
        "devpars",
        "subset",
        "plus",
        "ident",
        "kind"
    )

    fn = NULL
    default_devpars = NULL
    if ("ridge" %in% case$kind) {
        case$kind = "ridge"
        if (is.null(case$cols)) {
            case$cols = pal_ucscgb(alpha = .8)(26)
        }
        excluded_args = c(excluded_args, "split.by")
        fn = RidgePlot
        default_devpars = function(features, ncol, uidents) {
            if (is.null(ncol)) { ncol = 1 }
            list(
                width = 400 * ncol,
                height = ceiling(length(features) / ncol) * ifelse(length(uidents) < 10, 300, 400),
                res = 100
            )
        }
    } else if ("vln" %in% case$kind || "violin" %in% case$kind) {
        case$kind = "violin"
        if (is.null(case$cols)) {
            case$cols = pal_ucscgb(alpha = .8)(26)
        }
        fn = VlnPlot
        default_devpars = function(features, ncol, uidents) {
            if (is.null(ncol)) { ncol = 1 }
            list(
                width = 400 * ncol,
                height = ceiling(length(features) / ncol) * 400,
                res = 100
            )
        }
    } else if ("feature" %in% case$kind) {
        case$kind = "feature"
        if (is.null(case$cols)) {
            case$cols = c("lightgrey", pal_ucscgb()(1))
        }
        excluded_args = c(excluded_args, "group.by", "assay")
        case$shape.by = case$group.by
        fn = FeaturePlot
        default_devpars = function(features, ncol, uidents) {
            if (is.null(ncol)) { ncol = 1 }
            list(
                width = 400 * ncol,
                height = ceiling(length(features) / ncol) * 300,
                res = 100
            )
        }
    } else if ("dot" %in% case$kind) {
        case$kind = "dot"
        if (is.null(case$cols)) {
            case$cols = c("lightgrey", pal_ucscgb()(1))
        }
        if (is.null(case$plus)) {
            case$plus = 'theme_prism(axis_text_angle=90)'
        }
        excluded_args = c(excluded_args, "slot", "ncol")
        fn = DotPlot
        default_devpars = function(features, ncol, uidents) {
            list(
                height = max(length(uidents) * 80 + 150, 420),
                width = length(features) * 50 + 150,
                res = 100
            )
        }
    } else if ("heatmap" %in% case$kind) {
        case$kind = "heatmap"
        case = list_update(
            list(
                group.colors = pal_ucscgb(alpha = .8)(26),
                size = 3.5,
                group.bar.height = 0.01
            ),
            case
        )
        if (is.null(case$plus)) {
            case$plus = 'scale_fill_gradientn(colors = c("lightgrey", pal_ucscgb()(1)), na.value = "white")'
        }
        excluded_args = c(excluded_args, "group.by", "split.by", "downsample", "ncol")
        fn = DoHeatmap
        default_devpars = function(features, ncol, uidents) {
            list(
                width = length(uidents) * 60 + 150,
                height = length(features) * 40 + 150,
                res = 100
            )
        }
    } else if ("table" %in% case$kind) {
        case$kind = "table"
        excluded_args = c(excluded_args, "group.by", "split.by", "assay")
        case$assays = case$assay
        fn = AverageExpression
        if (is.null(case$slot)) {
            case$slot = "data"
        }
    } else {
        stop("Unknown kind of plot")
    }

    for (arg in excluded_args) {
        assign(arg, case[[arg]])
        case[[arg]] = NULL
    }

    if (is.character(subset)) {
        case$object = srtobj %>% filter(!!rlang::parse_expr(subset))
    } else {
        case$object = srtobj
    }
    if (!is.null(ident)) {
        Idents(case$object) = ident
    }
    case$features = .get_features(case$features)
    if (!is.null(case$ncol)) {
        case$ncol = min(case$ncol, length(case$features))
    }

    toc$kind = kind
    if (kind == "table") {
        expr = do_call(fn, case)$RNA %>%
            as.data.frame() %>%
            rownames_to_column("Feature") %>%
            select(Feature, everything())

        toc$file = paste0(slugify(name), ".txt")
        write.table(expr, file.path(odir, toc$file), sep="\t", quote=FALSE, row.names=FALSE)
    } else {
        devpars = list_update(
            default_devpars(case$features, case$ncol, unique(Idents(case$object))),
            devpars
        )
        if (kind == "heatmap") {
            if (!exists("downsample") || is.null(downsample)) {
                downsample = "average"
            }
            if (downsample %in% c("average", "mean")) {
                case$object = AverageExpression(case$object, return.seurat = TRUE)
            } else if (is.integer(downsample)) {
                case$object = subset(case$object, downsample = downsample)
            } else {
                stop("Unknown downsample method.")
            }
        }
        p = do_call(fn, case)
        if (!is.null(plus)) {
            for (pls in plus) {
                p = p + eval(parse(text = pls))
            }
        }
        figfile = file.path(odir, paste0(slugify(name), ".", kind, ".png"))
        toc$file = basename(figfile)
        png(figfile, width=devpars$width, height=devpars$height, res=devpars$res)
        tryCatch({
            print(p)
        }, error = function(e) {
            stop(
                paste(
                    paste(names(devpars), collapse=" "),
                    paste(devpars, collapse=" "),
                    e,
                    sep = "\n"
                )
            )
        })
        dev.off()
    }
    .add_toc(if (is.null(section)) name else section, toc)
}

sapply(names(features), do_one_features)
.save_toc()
