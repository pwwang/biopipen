# Loaded variables: srtfile, outdir, srtobj

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
    reduction <- if (reduction %in% c("dim", "auto")) scplotter:::default_dimreduc(srtobj) else reduction
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
    save_plot(p, prefix, devpars, selfcontained = FALSE)

    interactive_plot_file <- paste0(prefix, ".html")
    if (file.exists(interactive_plot_file)) {
        shared_asset_dir <- file.path(odir, "dimplot_3d_assets")
        asset_files_dir <- paste0(prefix, "_files")
        if (dir.exists(asset_files_dir)) {
            if (!dir.exists(shared_asset_dir)) {
                # copy the entire directory of assets to the shared asset directory
                dir.create(shared_asset_dir, recursive = TRUE, showWarnings = FALSE)
                file.copy(list.files(asset_files_dir, full.names = TRUE), shared_asset_dir, recursive = TRUE)
            }
            # remove the original directory of assets
            unlink(asset_files_dir, recursive = TRUE)

            # Now we need to update the paths in the HTML file to point to the shared asset directory
            html_content <- readLines(interactive_plot_file)
            updated_html_content <- gsub(paste0(basename(prefix), "_files/"), paste0("dimplot_3d_assets/"), html_content)
            writeLines(updated_html_content, interactive_plot_file)
        }
    }

    reporter$add(
        list(
            kind = "descr",
            content = paste0("Dimensionality reduction plot for ", case$group_by %||% "default identity")
        ),
        reporter$image(prefix, "pdf", FALSE),
        h1 = name
    )
}

invisible(sapply(names(dimplots), do_one_dimplot))
