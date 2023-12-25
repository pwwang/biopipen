# Loaded variables: srtfile, outdir, srtobj

ngenes_defaults <- {{envs.ngenes_defaults | r: todot="-"}}
ngenes <- {{envs.ngenes | r: todot="-", skip=1}}

odir <- file.path(outdir, "ngenes")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)

do_one_ngenes <- function(name) {
    log_info("Doing ngenes for: {name}")

    case <- list_update(ngenes_defaults, ngenes[[name]])
    case$devpars <- list_update(ngenes_defaults$devpars, case$devpars)

    figfile = file.path(odir, paste0(slugify(name), ".boxplot.png"))

    if (!is.null(case$subset)) {
        sobj <- srtobj %>% filter(!!rlang::parse_expr(case$subset))
    } else {
        sobj <- srtobj
    }
    df_cells <- sobj@meta.data %>% mutate(.nexpr = colSums(GetAssayData(sobj) > 0))

    select_cols = c(case$ident, case$group.by, case$split.by, ".nexpr")
    df_cells = df_cells %>% select(all_of(select_cols))

    p = df_cells %>%
        ggplot(aes(
            x=!!sym(case$ident),
            y=.nexpr,
            fill=!!sym(ifelse(is.null(case$group.by), case$ident, case$group.by))
        )) +
        geom_violin(position = ifelse(is.null(case$group.by), "identity", "dodge")) +
        geom_boxplot(
            position = ifelse(is.null(case$group.by), "identity", "dodge"),
            width = .1,
            fill = "white"
        ) +
        theme_prism(axis_text_angle = 90) +
        scale_fill_biopipen() +
        ylab("Number of genes expressed")

    if (!is.null(case$split.by)) {
        p = p + facet_wrap(case$split.by)
    }

    png(figfile, width=case$devpars$width, height=case$devpars$height, res=case$devpars$res)
    print(p)
    dev.off()

    add_report(
        list(
            kind = "descr",
            content = paste0(
                "Plots showing the number of genes expressed in each ",
                case$ident,
                ifelse(
                    is.null(case$group.by),
                    "",
                    paste0(", by ", paste0(case$group.by, collapse = ", "))
                )
            )
        ),
        list(kind = "image", src = figfile),
        h1 = name
    )
}

sapply(names(ngenes), do_one_ngenes)
