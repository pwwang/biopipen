# Loaded variables: srtfile, outdir, srtobj

dimplots_defaults = {{envs.dimplots_defaults | r: todot="-"}}
dimplots = {{envs.dimplots | r: todot="-", skip=1}}

odir = file.path(outdir, "dimplots")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)

do_one_dimplot = function(name) {
    log_info(paste0("Doing dimplots for: ", name))

    case = list_update(dimplots_defaults, dimplots[[name]])
    case$devpars = list_update(dimplots_defaults$devpars, dimplots[[name]]$devpars)
    if (!is.null(case$subset)) {
        case$object = srtobj %>% filter(!!rlang::parse_expr(case$subset))
    } else {
        case$object = srtobj
    }

    if (is.null(case$group.by)) { case$group.by = case$ident }
    key <- paste0("sub_umap_", case$ident)
    if (
        key %in% names(case$object@reductions) &&
        (is.null(case$reduction) || case$reduction %in% c("dim", "auto"))) {
        case$reduction = key
    }
    if (is.null(case$na_group)) {
        case$object = filter(case$object, !is.na(!!sym(case$group.by)))
    } else if (is.null(case$order)) {
        case$order = case$object@meta.data[[case$group.by]] %>%
            unique() %>% na.omit() %>% as.character() %>% sort()
        case$object@meta.data = replace_na(
            case$object@meta.data,
            setNames(list(case$na_group), case$group.by)
        )
    }
    case$cols = case$cols %||% pal_biopipen()(length(unique(case$object@meta.data[[case$group.by]])))
    if (!is.null(case$na_group)) {
        # Is the NA value in the first position?
        case$cols = c("lightgrey", case$cols[1:(length(case$cols) - 1)])
    }

    excluded_args = c("devpars", "ident", "subset")
    for (arg in excluded_args) {
        assign(arg, case[[arg]])
        case[[arg]] = NULL
    }

    if (case$reduction %in% c("dim", "auto")) { case$reduction = NULL }
    case$na_group <- NULL
    figfile = file.path(odir, paste0(slugify(name), ".dim.png"))
    png(figfile, width=devpars$width, height=devpars$height, res=devpars$res)
    p = do_call(DimPlot, case)
    print(p)
    dev.off()

    add_report(
        list(
            kind = "descr",
            content = paste0("Dimensionality reduction plot for ", case$group.by)
        ),
        list(
            kind = "image",
            src = figfile
        ),
        h1 = name
    )
}

sapply(names(dimplots), do_one_dimplot)
