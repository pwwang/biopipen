# Loaded variables: srtfile, outdir, srtobj

dimplots_defaults = {{envs.dimplots_defaults | r: todot="-"}}
dimplots = {{envs.dimplots | r: todot="-", skip=1}}

odir = file.path(outdir, "dimplots")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)

do_one_dimplot = function(name) {
    print(paste0("Doing dimplots for: ", name))

    case = list_update(dimplots_defaults, dimplots[[name]])
    case$devpars = list_update(dimplots_defaults$devpars, dimplots[[name]]$devpars)
    if (!is.null(case$subset)) {
        case$object = srtobj %>% filter(!!rlang::parse_expr(case$subset))
    } else {
        case$object = srtobj
    }

    if (is.null(case$group.by) && !is.null(case$ident)) {
        case$group.by = case$ident
    }

    n_uidents = length(unique(case$object@meta.data[[case$group.by]]))
    if (is.null(case$cols)) {
        case$cols = pal_biopipen()(n_uidents)
    }

    excluded_args = c("devpars", "ident", "subset")
    for (arg in excluded_args) {
        assign(arg, case[[arg]])
        case[[arg]] = NULL
    }

    if (case$reduction %in% c("dim", "auto")) {
        case$reduction = NULL
    }
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
