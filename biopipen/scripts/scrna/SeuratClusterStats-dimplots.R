# Loaded variables: srtfile, outdir, srtobj

# dimplots_defaults = {{envs.dimplots_defaults | r: todot="-"}}
# dimplots = {{envs.dimplots | r: todot="-", skip=1}}
log_info("dimplots:")

odir <- file.path(outdir, "dimplots")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)

do_one_dimplot = function(name) {
    log_info("- Case: {name}")

    case <- list_update(dimplots_defaults, dimplots[[name]])

    # Get functional arguments and inconsistent arguments
    use <- tolower(case$use); case$use <- NULL
    ident <- case$ident; case$ident <- NULL
    na.group <- case$na.group; case$na.group <- NULL
    show_na <- case$show_na; case$show_na <- NULL
    subset <- case$subset; case$subset <- NULL
    reduction <- case$reduction; case$reduction <- NULL
    label_repel <- case$label_repel; case$label_repel <- NULL
    repel <- case$repel; case$repel <- NULL
    devpars <- case$devpars; case$devpars <- NULL
    group.by <- case$group.by %||% ident; case$group.by <- NULL

    # Normalize arguments
    use <- match.arg(use, c("scp", "scp3d", "seurat"))
    reduction <- if (reduction %in% c("dim", "auto")) NULL else reduction
    devpars <- list_update(dimplots_defaults$devpars, devpars)
    key <- paste0("sub_umap_", ident)
    if (use == "scp" || use == "scp3d") {
        case$show_na <- show_na %||% !is.null(na.group)
        case$label_repel <- label_repel %||% repel
        case$group.by <- group.by
        if (!is.null(subset)) {
            case$srt <- srtobj %>% filter(!!parse_expr(subset))
        } else {
            case$srt <- srtobj
        }
        if (key %in% names(case$srt@reductions) && is.null(reduction)) {
            case$reduction = key
        }
        fun <- ifelse(use == "scp", CellDimPlot, CellDimPlot3D)
    } else {  # use == "seurat"
        case$repel <- repel %||% label_repel
        case$group.by <- group.by
        if (!is.null(subset)) {
            case$object <- srtobj %>% filter(!!parse_expr(subset))
        } else {
            case$object <- srtobj
        }
        na.group <- na.group %||% if (isTRUE(show_na)) "NA" else NULL
        if (is.null(na.group)) {
            case$object <- case$object %>% filter(!is.na(!!parse_expr(case$group.by)))
            case$cols <- case$cols %||% pal_biopipen()(length(unique(case$object@meta.data[[case$group.by]])))
        } else {
            case$order = case$object@meta.data[[case$group.by]] %>%
                unique() %>% na.omit() %>% as.character() %>% sort()
            case$object@meta.data = replace_na(
                case$object@meta.data,
                setNames(list(na_group), case$group.by)
            )
            # Is the NA value in the first position?
            case$cols <- c("lightgrey", case$cols[1:(length(case$cols) - 1)])
        }
        if (key %in% names(case$object@reductions) && is.null(reduction)) {
            case$reduction = key
        }
        fun <- DimPlot
    }
    case <- case[names(case) %in% formalArgs(fun)]

    if (use == "scp3d") {
        case$save <- file.path(odir, paste0(slugify(name), ".dim.html"))
        do_call(fun, case)

        add_report(
            list(
                kind = "descr",
                content = paste0("Dimensionality reduction plot for ", case$group.by)
            ),
            list(
                kind = "tag",
                tag = "Plotly",
                title = paste0("Dimensionality reduction plot for ", case$group.by),
                width = devpars$width,
                height = devpars$height,
                src = case$save
            ),
            h1 = name
        )
    } else {
        figfile <- file.path(odir, paste0(slugify(name), ".dim.png"))
        png(figfile, width=devpars$width, height=devpars$height, res=devpars$res)
        p <- do_call(fun, case)
        print(p)
        dev.off()
        pdffile <- file.path(odir, paste0(slugify(name), ".dim.pdf"))
        ggsave(pdffile, p, width=devpars$width, height=devpars$height, dpi=devpars$res, device="pdf", units="px")

        add_report(
            list(
                kind = "descr",
                content = paste0("Dimensionality reduction plot for ", case$group.by)
            ),
            list(
                kind = "image",
                src = figfile,
                download = pdffile
            ),
            h1 = name
        )
    }
}

sapply(names(dimplots), do_one_dimplot)
