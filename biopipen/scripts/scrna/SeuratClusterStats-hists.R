# Loaded variables: srtfile, outdir, srtobj

hists_defaults <- {{envs.hists_defaults | r: todot="-"}}
hists <- {{envs.hists | r: todot="-", skip=1}}

do_one_hists <- function(m, case, odir, h1, each = NULL) {
    ofile <- file.path(odir, paste0(slugify(h1), ifelse(is.null(each), "", paste0("-", slugify(each))), ".png"))

    p <- ggplot(m, aes(x=!!sym(case$x))) +
        geom_histogram(
            bins = case$bins,
            fill = "white",
            color = "#3a9ddc",
            stat = ifelse(is.numeric(m[[case$x]]), "bin", "count")
        ) +
        theme_prism() +
        xlab(case$x) +
        ylab("Number of cells")

    if (!is.null(case$cells_by)) {
        p <- p + facet_wrap(vars(!!sym(case$cells_by)), ncol = case$ncol, drop = FALSE)
        ngroups <- length(unique(m[[case$cells_by]]))
    } else {
        ngroups <- 1
    }

    if (!is.null(case$plus) && length(case$plus) > 0) {
        for (pl in case$plus) {
            p <- p + eval(parse(text=pl))
        }
    }

    devpars <- case$devpars
    if (is.null(devpars$res)) { devpars$res <- 100 }
    if (is.null(devpars$width)) {
        devpars$width <- ifelse(ngroups == 1, 600, 300 * min(ngroups, case$ncol)) }
    if (is.null(devpars$height)) {
        devpars$height <- ifelse(ngroups == 1, 300, 150 * ceiling(ngroups / case$ncol))
    }
    png(ofile, width=devpars$width, height=devpars$height, res=devpars$res)
    print(p)
    dev.off()

    # Add report
    if (!is.null(each)) {
        add_report(
            list(src = ofile, descr = first(m[[case$each]])),
            h1 = h1,
            ui = "table_of_images"
        )
    } else {
        add_report(
            list(kind = "image", src = ofile),
            h1 = h1
        )
    }
}

if (is.null(hists) || length(hists) == 0) {
    log_warn("No hists cases specified, skipping ...")
} else {

    for (name in names(hists)) {
        hists[[name]] <- list_update(hists_defaults, hists[[name]])
        case <- hists[[name]]

        odir <- file.path(outdir, "hists", slugify(name))
        dir.create(odir, recursive=TRUE, showWarnings=FALSE)

        h1 <- name
        if (!is.null(case$subset)) {
            meta <- srtobj@meta.data %>% filter(!!rlang::parse_expr(case$subset))
        } else {
            meta <- srtobj@meta.data
        }

        if (is.null(case$x)) { stop("`x` must be specified for hists.") }

        meta <- meta %>% drop_na(!!sym(case$x))
        if (!is.null(case$cells_by)) { meta <- meta %>% drop_na(!!sym(case$cells_by)) }
        if (!is.null(case$each)) { meta <- meta %>% drop_na(!!sym(case$each)) }
        if (!is.null(case$x_order) && length(case$x_order) > 0) {
            meta <- meta %>% filter(!!sym(case$x) %in% case$x_order) %>%
                arrange(match(!!sym(case$x), case$x_order))
        }
        if (!is.null(case$cells_order) && length(case$cells_order) > 0) {
            if (is.null(case$cells_by)) { stop("`cells_by` must be specified for `cells_order`.") }
            meta <- meta %>% filter(!!sym(case$cells_by) %in% case$cells_order) %>%
                arrange(match(!!sym(case$cells_by), case$cells_order))
        } else if (!is.null(case$cells_orderby)) {
            meta <- meta %>% arrange(!!parse_expr(case$cells_orderby))
            cell_groups <- meta %>% pull(!!sym(case$cells_by)) %>% unique()
            meta <- meta %>% filter(!!sym(case$cells_by) %in% cell_groups[1:case$cells_n])
            meta[[case$cells_by]] <- factor(meta[[case$cells_by]], levels=cell_groups[1:case$cells_n])
        }

        if (!is.null(case$each)) {
            eachs <- meta %>% pull(!!sym(case$each)) %>% unique()
            add_report(
                list(
                    kind = "descr",
                    content = paste0(
                        "Histogram showing number of cells ",
                        ifelse(is.null(case$cells_by), "", paste0("by <code>", case$cells_by, "</code> ")),
                        "in each <code>",
                        case$each,
                        "</code> along <code>",
                        case$x,
                        "</code>"
                    )
                ),
                h1 = h1
            )
            for (each in eachs) {
                log_info("Doing hists for: {h1} - {each} ...")
                m <- meta %>% filter(!!sym(case$each) == each)
                do_one_hists(m, case, odir, h1, each)
            }
        } else {
            log_info("Doing hists for: {h1} ...")
            add_report(
                list(
                    kind = "descr",
                    content = paste0(
                        "Histogram showing number of cells ",
                        ifelse(is.null(case$cells_by), "", paste0("by <code>", case$cells_by, "</code> ")),
                        "along <code>",
                        case$x,
                        "</code>"
                    )
                ),
                h1 = h1
            )
            do_one_hists(meta, case, odir, h1)
        }
    }

}