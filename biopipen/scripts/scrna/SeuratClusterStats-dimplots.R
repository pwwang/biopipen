# Loaded variables: srtfile, outdir, srtobj

# dimplots_defaults = {{envs.dimplots_defaults | r: todot="-"}}
# dimplots = {{envs.dimplots | r: todot="-", skip=1}}
log$info("dimplots:")

odir <- file.path(outdir, "dimplots")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)

do_one_dimplot = function(name) {
    log$info("- Case: {name}")

    case <- list_update(dimplots_defaults, dimplots[[name]])

    # Get functional arguments and inconsistent arguments
    subset <- case$subset; case$subset <- NULL
    reduction <- case$reduction; case$reduction <- NULL
    devpars <- case$devpars; case$devpars <- NULL
    stopifnot("Either 'group_by' or 'ident' should be specified in dimplots, not both." =
        is.null(case$group_by) || is.null(case$ident) || identical(case$group_by, case$ident))

    # Normalize arguments
    reduction <- if (reduction %in% c("dim", "auto")) DefaultDimReduc(srtobj) else reduction
    devpars <- list_update(dimplots_defaults$devpars, devpars)

    case$group_by <- case$group_by %||% case$ident %||% GetIdentityColumn(srtobj)
    # key <- paste0("sub_umap_", case$group_by)
    subcluster_key <- paste0(case$group_by, ".", reduction)

    if (!is.null(subset)) {
        case$object <- srtobj %>% filter(!!parse_expr(subset))
    } else {
        case$object <- srtobj
    }
    if (subcluster_key %in% names(case$object@reductions)) {
        case$reduction = subcluster_key
    } else {
        case$reduction = reduction
    }

    p <- do_call(CellDimPlot, case)
    prefix <- file.path(odir, paste0(slugify(name), ".dim"))
    save_plot(p, prefix, devpars)

    reporter$add(
        list(
            kind = "descr",
            content = paste0("Dimensionality reduction plot for ", case$group_by)
        ),
        reporter$image(prefix, "pdf", FALSE),
        h1 = name
    )
}

sapply(names(dimplots), do_one_dimplot)
