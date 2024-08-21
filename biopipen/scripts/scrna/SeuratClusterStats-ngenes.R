# Loaded variables: srtfile, outdir, srtobj

# ngenes_defaults <- {{envs.ngenes_defaults | r: todot="-"}}
# ngenes <- {{envs.ngenes | r: todot="-", skip=1}}
log_info("ngenes:")

odir <- file.path(outdir, "ngenes")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)

do_one_ngenes <- function(name) {
    log_info("- Case: {name}")

    case <- list_update(ngenes_defaults, ngenes[[name]])
    case$devpars <- list_update(ngenes_defaults$devpars, case$devpars)

    if (!is.null(case$subset)) {
        sobj <- srtobj %>% filter(!!rlang::parse_expr(case$subset))
    } else {
        sobj <- srtobj
    }
    df_cells <- sobj@meta.data %>% mutate(.nexpr = colSums(GetAssayData(sobj) > 0))

    select_cols = c(case$ident, case$group.by, case$split.by, ".nexpr")
    df_cells = df_cells %>% select(all_of(select_cols))

    p = df_cells |>
        ggplot(aes(
            x=!!sym(case$ident),
            y=.nexpr,
            fill=!!sym(ifelse(is.null(case$group.by), case$ident, case$group.by))
        )) +
        geom_violin(alpha = 0.6, position = ifelse(is.null(case$group.by), "identity", "dodge")) +
        geom_boxplot(
            position = ifelse(is.null(case$group.by), "identity", "dodge"),
            width = .1,
            fill = "white"
        ) +
        theme_prism(axis_text_angle = 45) +
        scale_fill_biopipen() +
        ylab("Number of genes expressed")

    if (!is.null(case$split.by)) {
        p = p + facet_wrap(case$split.by)
    }

    figprefix = file.path(odir, paste0(slugify(name), ".boxplot"))

    save_plot(p, figprefix, case$devpars)
    save_plotcode(
        p,
        c(
            'library(rlang)',
            'library(ggplot2)',
            'library(ggprism)',
            '',
            'load("data.RData")'
        ),
        figprefix,
        "df_cells", "case", "scale_fill_biopipen", "pal_biopipen"
    )

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
        list(
            kind = "image",
            src = paste0(figprefix, ".png"),
            download = list(
                paste0(figprefix, ".pdf"),
                list(
                    src = paste0(figprefix, ".code.zip"),
                    tip = "Download the code to reproduce the plot",
                    icon = "Code"
                )
            )),
        h1 = name
    )
}

sapply(names(ngenes), do_one_ngenes)
