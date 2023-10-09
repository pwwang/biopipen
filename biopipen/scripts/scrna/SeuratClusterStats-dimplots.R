# Loaded variables: srtfile, outdir, srtobj

dimplots_defaults = {{envs.dimplots_defaults | r: todot="-"}}
dimplots = {{envs.dimplots | r: todot="-", skip=1}}

odir = file.path(outdir, "dimplots")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)
report_toc_file = file.path(odir, "report_toc.json")
# Realname => file
report_toc = list()

do_one_dimplot = function(name) {
    print(paste0("Doing dimplots for: ", name))

    case = list_update(dimplots_defaults, dimplots[[name]])
    case$devpars = list_update(dimplots_defaults$devpars, dimplots[[name]]$devpars)
    case$object = srtobj
    if (is.null(case$cols)) {
        case$cols = pal_ucscgb()(26)
    }

    excluded_args = c("devpars", "ident")
    for (arg in excluded_args) {
        assign(arg, case[[arg]])
        case[[arg]] = NULL
    }

    if (case$reduction %in% c("dim", "auto")) {
        case$reduction = NULL
    }
    report_toc[[name]] <<- paste0(slugify(name), ".dim.png")
    figfile = file.path(odir, report_toc[[name]])
    png(figfile, width=devpars$width, height=devpars$height, res=devpars$res)
    p = do_call(DimPlot, case)
    print(p)
    dev.off()
}

sapply(names(dimplots), do_one_dimplot)
.save_toc()
