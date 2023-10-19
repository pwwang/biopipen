# Loaded variables: srtfile, outdir, srtobj

stats_defaults = {{envs.stats_defaults | r: todot="-"}}
stats = {{envs.stats | r: todot="-", skip=1}}

odir = file.path(outdir, "stats")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)
report_toc_file = file.path(odir, "report_toc.json")
# Realname => {bar: ..., pie: ..., table: ...}
report_toc = list()

.add_toc = function(name, toc) {
    report_toc[[name]] <<- toc
}

.save_toc = function() {
    writeLines(toJSON(report_toc, pretty = TRUE, auto_unbox = TRUE), report_toc_file)
}

do_one_stats = function(name) {
    print(paste0("Doing stats for: ", name))

    toc = list()

    case = list_update(stats_defaults, stats[[name]])
    case$devpars = list_update(stats_defaults$devpars, case$devpars)
    if (isTRUE(case$pie) && !is.null(case$group.by)) {
        stop(paste0(name, ": pie charts are not supported for group-by"))
    }

    figfile = file.path(odir, paste0(slugify(name), ".bar.png"))
    piefile = file.path(odir, paste0(slugify(name), ".pie.png"))
    tablefile = file.path(odir, paste0(slugify(name), ".txt"))

    df_cells = srtobj@meta.data
    if (!is.null(case$subset)) {
        df_cells = df_cells %>% filter(!!rlang::parse_expr(case$subset))
    }

    select_cols = c(case$ident, case$group.by, case$split.by)
    df_cells = df_cells %>%
        select(all_of(select_cols)) %>%
        group_by(!!!syms(select_cols)) %>%
        summarise(.n = n(), .groups = "drop")  %>%
        mutate(.frac = .n / sum(.n))

    if (isTRUE(case$table)) {
        toc$table = basename(tablefile)
        write.table(df_cells, tablefile, sep="\t", quote=FALSE, row.names=FALSE)
    }
    if (isTRUE(case$pie)) {
        p_pie = df_cells %>%
            arrange(!!sym(case$ident)) %>%
            ggplot(aes(x="", y=.n, fill=!!sym(case$ident))) +
            geom_bar(stat="identity", width=1, alpha=.8, position = position_stack(reverse = TRUE)) +
            coord_polar("y", start=0) +
            scale_fill_ucscgb(alpha=.8) +
            guides(fill = guide_legend(title = case$ident)) +
            theme_void() +
            geom_label(
                if (isTRUE(case$frac))
                    aes(label=sprintf("%.1f%%", .frac * 100))
                else
                    aes(label=.n),
                position = position_stack(vjust = 0.5),
                color="#333333",
                fill="#EEEEEE",
                size=5
            )

        if (!is.null(case$split.by)) {
            p_pie = p_pie + facet_wrap(case$split.by)
        }

        toc$pie = basename(piefile)
        png(piefile, width=case$devpars$width, height=case$devpars$height, res=case$devpars$res)
        print(p_pie)
        dev.off()
    }

    ngroups = ifelse(is.null(case$group.by), 1, length(unique(df_cells[[case$group.by]])))
    nidents = length(unique(df_cells[[case$ident]]))
    bar_position = ifelse(ngroups > 5, "stack", "dodge")
    p = df_cells %>%
        ggplot(aes(
            x=!!sym(case$ident),
            y=if (isTRUE(case$frac)) .frac else .n,
            fill=!!sym(ifelse(is.null(case$group.by), case$ident, case$group.by))
        )) +
        geom_bar(stat="identity", position=bar_position, alpha=.8) +
        theme_prism(axis_text_angle = 90) +
        scale_fill_manual(values=rep(pal_ucscgb(alpha=.8)(26), 10)[1:max(ngroups, nidents)]) +
        ylab(ifelse(isTRUE(case$frac), "Fraction of cells", "Number of cells"))

    if (!is.null(case$split.by)) {
        p = p + facet_wrap(case$split.by)
    }

    toc$bar = basename(figfile)
    png(figfile, width=case$devpars$width, height=case$devpars$height, res=case$devpars$res)
    print(p)
    dev.off()

    .add_toc(name, toc)
}

sapply(names(stats), do_one_stats)
.save_toc()
