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
        "cluster_orderby",
        "kind"
    )

    if (is.character(case$subset)) {
        case$object = srtobj %>% filter(!!rlang::parse_expr(case$subset))
    } else {
        case$object = srtobj
    }

    if (!is.null(case$ident)) {
        case$object = case$object %>% filter(!is.na(!!sym(case$ident)))
        Idents(case$object) = case$ident
    }
    cluster_order_val <- NULL
    if (!is.null(case$cluster_orderby) && length(case$cluster_orderby) > 1) {
        case$object@meta.data[[case$ident]] = factor(Idents(case$object), levels = case$cluster_orderby)
        Idents(case$object) = case$ident
    } else if (!is.null(case$cluster_orderby)) {
        cluster_order_df = case$object@meta.data %>%
            group_by(!!sym(case$ident)) %>%
            summarise(!!sym(case$cluster_orderby) := !!parse_expr(case$cluster_orderby)) %>%
            arrange(!!sym(case$cluster_orderby))
        cluster_order_val = cluster_order_df[[case$cluster_orderby]]
        clusters = cluster_order_df[[case$ident]] %>% as.character() %>% unique()
        case$object@meta.data[[case$ident]] = factor(Idents(case$object), levels = clusters)
        Idents(case$object) = case$ident
    }
    n_uidents = length(unique(Idents(case$object)))
    max_nchar_idents = max(nchar(unique(as.character(Idents(case$object)))))

    fn = NULL
    default_devpars = NULL
    if (case$kind %in% c("ridge", "ridgeplot")) {
        case$kind = "ridge"
        if (is.null(case$cols)) {
            case$cols = pal_biopipen()(n_uidents)
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
        if (is.null(case$cols)) { case$cols = pal_biopipen()(n_uidents) }
        if (is.null(case$pt.size)) { case$pt.size = 0 }

        excluded_args = c(excluded_args, "reduction")
        fn = VlnPlot
        default_devpars = function(features, ncol) {
            if (is.null(ncol)) { ncol = 1 }
            list(
                width = 400 * ncol,
                height = ceiling(length(features) / ncol + max_nchar_idents * .05) * 150,
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
            case$plus = 'scale_fill_gradient2(
                low = "lightblue",
                high = "darkblue",
                na.value = "white"
            )'
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
    } else if (case$kind == "avgheatmap") {
        case$kind = "avgheatmap"
        excluded_args = c(
            excluded_args,
            "group.by", "split.by", "downsample", "ncol", "reduction", "layer",
            "assay", "object"
        )
        default_devpars = function(features, ncol) {
            list(
                width = n_uidents * 30 + 350,
                height = length(features) * 15 + 150,
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
    if (kind == "avgheatmap") {
        figfile <- file.path(odir, paste0(slugify(name), ".avgheatmap.png"))
        assay <- assay %||% DefaultAssay(object)
        layer <- layer %||% ifelse("scale.data" %in% Layers(object, assay = assay), "scale.data", "data")

        case_features <- case$features
        case$features <- NULL
        meta_feats <- intersect(case_features, colnames(object@meta.data))
        expr_feats <- setdiff(case_features, meta_feats)
        exprs <- NULL
        if (length(meta_feats) > 0) {
            exprs <- object@meta.data %>% select(all_of(c(meta_feats, ident))) %>%
                group_by(!!sym(ident)) %>%
                summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
                column_to_rownames(ident) %>%
                t()
        }
        if (length(expr_feats)) {
            exprs_tmp <- AverageExpression(
                object,
                assays = assay,
                layer = layer,
                features = expr_feats,
                group.by = ident)[[assay]]
            exprs <- bind_rows(exprs, as.data.frame(exprs_tmp))
        }

        ha <- NULL
        extra_height <- 0
        extra_width <- 0  # legend
        if (!is.null(cluster_order_val)) {
            ha <- list()
            ha[[cluster_orderby]] <- cluster_order_val
            if (is.numeric(cluster_order_val)) {
                col_fun <- colorRamp2(
                    c(min(cluster_order_val), max(cluster_order_val)),
                    c("lightyellow", "red"))
                ha$col <- list()
                ha$col[[cluster_orderby]] <- col_fun
            }
            ha <- do_call(HeatmapAnnotation, ha)
            extra_height <- 40
            extra_width <- 120
        }

        col_fun <- colorRamp2(
            c(min(exprs, na.rm = T), 0, max(exprs, na.rm = T)),
            c("lightblue", "white", "darkred"))

        case <- list_update(list(
            matrix = as.matrix(exprs),
            name = "Average expression",
            col = col_fun,
            na_col = "white",
            row_names_side = "right",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            rect_gp = gpar(col = "gray", lwd = 1),
            row_names_max_width = max_text_width(rownames(exprs)),
            top_annotation = ha
        ), case)
        p <- do_call(Heatmap, case)

        def_devpars = default_devpars(case_features, NULL)
        def_devpars$width = def_devpars$width + extra_width
        def_devpars$height = def_devpars$height + extra_height
        devpars = list_update(def_devpars, devpars)
        png(figfile, res = devpars$res, width = devpars$width, height = devpars$height)
        print(p)
        dev.off()

        add_report(
            list(
                kind = "descr",
                content = paste0("Average expression values for selected features, by ", ident)
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
            log_warn("- `downsample` is not specified for `heatmap`, using `downsample=1000`")
            downsample = 1000
        }
        if (is.numeric(downsample)) {
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
