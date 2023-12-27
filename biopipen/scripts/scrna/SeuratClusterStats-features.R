# Loaded variables: srtfile, outdir, srtobj

features_defaults = {{envs.features_defaults | r: todot="-"}}
features = {{envs.features | r: todot="-", skip=1}}

odir = file.path(outdir, "features")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)

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
    log_info("Doing features for: {name}")

    case = list_update(features_defaults, features[[name]])
    case$devpars = list_update(features_defaults$devpars, features[[name]]$devpars)
    excluded_args = c(
        "section",
        "devpars",
        "subset",
        "plus",
        "ident",
        "kind"
    )

    if (is.character(case$subset)) {
        case$object = srtobj %>% filter(!!rlang::parse_expr(case$subset))
    } else {
        case$object = srtobj
    }

    if (!is.null(case$ident)) {
        Idents(case$object) = case$ident
    }
    n_uidents = length(unique(Idents(case$object)))

    fn = NULL
    default_devpars = NULL
    if (case$kind %in% c("ridge", "ridgeplot")) {
        case$kind = "ridge"
        if (is.null(case$cols)) {
            case$cols = pal_biopipen()(32)
        }
        excluded_args = c(excluded_args, "split.by", "reduction")
        fn = RidgePlot
        default_devpars = function(features, ncol) {
            if (is.null(ncol)) { ncol = 1 }
            list(
                width = 400 * ncol,
                height = ceiling(length(features) / ncol) * ifelse(n_uidents < 10, 300, 400),
                res = 100
            )
        }
    } else if (case$kind %in% c("vln", "violin", "vlnplot", "violinplot")) {
        case$kind = "violin"
        if (is.null(case$cols)) {
            case$cols = pal_biopipen()(n_uidents)
        }
        excluded_args = c(excluded_args, "reduction")
        fn = VlnPlot
        default_devpars = function(features, ncol) {
            if (is.null(ncol)) { ncol = 1 }
            list(
                width = 400 * ncol,
                height = ceiling(length(features) / ncol) * 400,
                res = 100
            )
        }
    } else if (case$kind %in% c("feature", "featureplot")) {
        case$kind = "feature"
        if (is.null(case$cols)) {
            case$cols = c("lightgrey", pal_biopipen()(1))
        }
        excluded_args = c(excluded_args, "group.by", "assay", "layer")
        case$shape.by = case$group.by
        if (!is.null(case$ident)) {
            key <- paste0("sub_umap_", case$ident)
            if (key %in% names(case$object@reductions) && is.null(case$reduction)) {
                case$reduction = key
                case$object = filter(case$object, !is.na(!!sym(case$ident)))
            }
        }
        fn = FeaturePlot
        default_devpars = function(features, ncol) {
            if (is.null(ncol)) { ncol = 1 }
            list(
                width = 400 * ncol,
                height = ceiling(length(features) / ncol) * 300,
                res = 100
            )
        }
    } else if (case$kind %in% c("dot", "dotplot")) {
        case$kind = "dot"
        if (is.null(case$cols)) {
            case$cols = c("lightgrey", pal_biopipen()(1))
        }
        if (is.null(case$plus)) {
            case$plus = 'theme_prism(axis_text_angle=90)'
        }
        excluded_args = c(excluded_args, "layer", "ncol", "reduction")
        fn = DotPlot
        default_devpars = function(features, ncol) {
            list(
                height = max(n_uidents * 80 + 150, 420),
                width = length(features) * 50 + 150,
                res = 100
            )
        }
    } else if ("heatmap" == case$kind) {
        case$kind = "heatmap"
        case = list_update(
            list(
                group.colors = pal_biopipen()(n_uidents),
                size = 3.5,
                group.bar.height = 0.01
            ),
            case
        )
        if (is.null(case$plus)) {
            case$plus = 'scale_fill_gradientn(colors = c("lightgrey", pal_biopipen()(1)), na.value = "white")'
        }
        excluded_args = c(excluded_args, "group.by", "split.by", "downsample", "ncol", "reduction", "layer")
        fn = DoHeatmap
        default_devpars = function(features, ncol) {
            list(
                width = n_uidents * 60 + 150,
                height = length(features) * 40 + 150,
                res = 100
            )
        }
    } else if (case$kind %in% c("bar", "barplot")) {
        case$kind = "bar"
        if (is.null(case$features) || length(case$features) == 0) {
            stop("No features is specified for barplot")
        }
        if (length(case$features) > 1) {
            stop("Only one feature is allowed for barplot")
        }
        excluded_args = c(excluded_args, "reduction")
        default_devpars = function(features, ncol) {
            if (is.null(ncol)) { ncol = 1 }
            list(
                height = 500 * ncol,
                width = n_uidents * 60 + 150,
                res = 100
            )
        }
    } else if ("table" == case$kind) {
        case$kind = "table"
        excluded_args = c(excluded_args, "group.by", "split.by", "assay", "reduction")
        case$assays = case$assay
        fn = AverageExpression
        if (is.null(case$layer)) { case$layer = "data" }
    } else {
        stop(paste0("Unknown kind of plot: ", case$kind))
    }

    for (arg in excluded_args) {
        assign(arg, case[[arg]])
        case[[arg]] = NULL
    }

    if (kind == "bar") {
        figfile <- file.path(odir, paste0(slugify(name), ".bar.png"))
        genes <- rownames(GetAssayData(case$object))
        genes <- genes[sapply(genes, function(x) grepl(x, case$features))]
        if (length(genes) == 0) {
            p <- case$object@meta.data %>%
                group_by(Idents = Idents(case$object)) %>%
                summarise(!!sym(name) := !!parse_expr(case$features)) %>%
                ggplot(aes(x = Idents, y = !!sym(name)))
        } else {
            p <- case$object@meta.data %>%
                bind_cols(FetchData(case$object, vars = genes, layer = case$layer, cells = rownames(case$object@meta.data))) %>%
                group_by(Idents = Idents(case$object)) %>%
                summarise(!!sym(name) := !!parse_expr(case$features)) %>%
                ggplot(aes(x = Idents, y = !!sym(name)))
        }

        if (!is.null(case$group.by)) {
            p <- p + geom_bar(aes(fill = !!sym(case$group.by)), stat = "identity", position = "dodge")
        } else {
            p <- p + geom_bar(aes(fill = Idents), stat = "identity")
        }

        p <- p +
            scale_fill_biopipen() +
            theme_prism() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            labs(x = "Idents", y = name)

        if (!is.null(case$split.by)) {
            p <- p + facet_wrap(~ !!sym(case$split.by), ncol = case$ncol)
        } else {
            case$ncol = 1
        }

        if (!is.null(case$plus)) {
            p <- p + eval(parse(text = case$plus))
        }
        devpars = list_update(default_devpars(NULL, case$ncol), devpars)
        png(figfile, res = devpars$res, width = devpars$width, height = devpars$height)
        print(p)
        dev.off()

        add_report(
            list(
                kind = "descr",
                content = paste0(kind, "plots for selected features <code>", case$features, "</code>, by ", ident)
            ),
            list(
                kind = "image",
                src = figfile
            ),
            h1 = ifelse(is.null(section), name, section),
            h2 = ifelse(is.null(section), "#", name)
        )

        return(NULL)
    }

    case$features = .get_features(case$features)
    if (!is.null(case$ncol)) {
        case$ncol = min(case$ncol, length(case$features))
    }

    if (kind == "table") {
        expr = do_call(fn, case)[[DefaultAssay(case$object)]] %>%
            as.data.frame() %>%
            rownames_to_column("Feature") %>%
            select(Feature, everything())

        exprfile = file.path(odir, paste0(slugify(name), ".txt"))
        write.table(expr, exprfile, sep="\t", quote=FALSE, row.names=FALSE)

        add_report(
            list(
                kind = "descr",
                content = paste0("Table of expression value for selected features, by ", ident)
            ),
            list(
                kind = "table",
                src = exprfile
            ),
            h1 = ifelse(is.null(case$section), name, case$section),
            h2 = ifelse(is.null(case$section), "#", name)
        )
        return(NULL)
    }

    devpars = list_update(default_devpars(case$features, case$ncol), devpars)
    if (kind == "heatmap") {
        if (!exists("downsample") || is.null(downsample)) {
            downsample = "average"
        }
        if (downsample %in% c("average", "mean")) {
            case$object = AverageExpression(case$object, return.seurat = TRUE)
        } else if (is.integer(downsample)) {
            case$object = base::subset(case$object, downsample = downsample)
        } else {
            stop(paste0("Unknown downsample method: ", downsample))
        }
    }
    p = do_call(fn, case)
    if (!is.null(plus)) {
        for (pls in plus) {
            p = p + eval(parse(text = pls))
        }
    }
    figfile = file.path(odir, paste0(slugify(name), ".", slugify(case$kind), ".png"))
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

    add_report(
        list(
            kind = "descr",
            content = paste0(kind, "plots for selected features, by ", ident)
        ),
        list(
            kind = "image",
            src = figfile
        ),
        h1 = ifelse(is.null(section), name, section),
        h2 = ifelse(is.null(section), "#", name)
    )
}

sapply(names(features), do_one_features)
