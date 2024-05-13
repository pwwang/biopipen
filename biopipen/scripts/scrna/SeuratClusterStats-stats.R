# Loaded variables: srtfile, outdir, srtobj
library(circlize)

stats_defaults = {{envs.stats_defaults | r: todot="-"}}
stats = {{envs.stats | r: todot="-", skip=1}}

odir = file.path(outdir, "stats")
dir.create(odir, recursive=TRUE, showWarnings=FALSE)

do_one_stats = function(name) {
    log_info("Doing stats for: {name}")

    case = list_update(stats_defaults, stats[[name]])
    case$devpars = list_update(stats_defaults$devpars, case$devpars)
    case$pie_devpars = list_update(stats_defaults$pie_devpars, case$pie_devpars)
    case$circos_devpars = list_update(stats_defaults$circos_devpars, case$circos_devpars)
    if (isTRUE(case$pie) && !is.null(case$group.by)) {
        stop(paste0(name, ": pie charts are not supported for group-by"))
    }
    if (!isTRUE(case$frac) && isTRUE(case$frac_ofall)) {
        stop(paste0(name, ": frac_ofall is only supported when frac is true"))
    }
    if (isTRUE(case$frac_ofall) && is.null(case$group.by)) {
        stop(paste0(name, ": frac_ofall is only supported for group-by"))
    }
    if (isTRUE(case$transpose) && is.null(case$group.by)) {
        stop(paste0(name, ": transpose is only supported for group-by"))
    }

    figfile = file.path(odir, paste0(slugify(name), ".bar.png"))
    piefile = file.path(odir, paste0(slugify(name), ".pie.png"))
    circosfile = file.path(odir, paste0(slugify(name), ".circos.png"))
    samtablefile = file.path(odir, paste0(slugify(name), ".bysample.txt"))
    tablefile = file.path(odir, paste0(slugify(name), ".txt"))

    df_cells = srtobj@meta.data %>% drop_na(!!sym(case$ident))
    if (!is.null(case$subset)) {
        df_cells = df_cells %>% filter(!!rlang::parse_expr(case$subset))
    }

    select_cols = unique(c(case$ident, case$group.by, case$split.by))
    if (!is.null(case$split.by)) {
        plot_df = do_call(rbind, lapply(group_split(
            df_cells %>% select(all_of(select_cols)),
            !!!syms(case$split.by)
        ), function(df) {
            out <- df %>% group_by(!!!syms(select_cols)) %>% summarise(.n = n(), .groups = "drop")
            if (!is.null(case$group.by) && isTRUE(case$frac)) {
                if (isTRUE(case$frac_ofall)) {
                    out <- out %>% mutate(.frac = .n / sum(.n))
                } else if (isTRUE(case$transpose)) {
                    out <- out %>% group_by(!!sym(case$ident)) %>% mutate(.frac = .n / sum(.n))
                } else {
                    out <- out %>% group_by(!!sym(case$group.by)) %>% mutate(.frac = .n / sum(.n))
                }
            }
            out
        }))
    } else if (!is.null(case$group.by) && isTRUE(case$frac)) {
        plot_df <- df_cells %>%
            select(all_of(select_cols)) %>%
            group_by(!!!syms(select_cols)) %>%
            summarise(.n = n(), .groups = "drop")
        if (isTRUE(case$frac_ofall)) {
            plot_df = plot_df %>% mutate(.frac = .n / sum(.n))
        } else {
            plot_df = plot_df %>%
                group_by(!!sym(ifelse(isTRUE(case$transpose), case$group.by, case$ident))) %>%
                mutate(.frac = .n / sum(.n))
        }
    } else {
        plot_df <- df_cells %>%
            select(all_of(select_cols)) %>%
            group_by(!!!syms(select_cols)) %>%
            summarise(.n = n(), .groups = "drop")

        if (isTRUE(case$frac) || isTRUE(case$frac_ofall)) {
            plot_df <- plot_df %>% mutate(.frac = .n / sum(.n))
        }
    }

    write.table(plot_df, tablefile, sep="\t", quote=FALSE, row.names=FALSE)

    ngroups = ifelse(is.null(case$group.by), 1, length(unique(plot_df[[case$group.by]])))
    nidents = length(unique(plot_df[[case$ident]]))
    bar_position = ifelse(case$position == "auto", ifelse(ngroups > 5, "stack", "dodge"), case$position)
    p = plot_df %>%
        ggplot(aes(
            x=!!sym(ifelse(case$transpose, case$group.by, case$ident)),
            y=if (isTRUE(case$frac)) .frac else .n,
            fill=!!sym(ifelse(is.null(case$group.by) || isTRUE(case$transpose), case$ident, case$group.by))
        )) +
        geom_bar(stat="identity", position=bar_position, alpha=.8) +
        theme_prism(axis_text_angle = 90) +
        scale_fill_biopipen() +
        ylab(ifelse(isTRUE(case$frac), "Fraction of cells", "Number of cells"))

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
                "Plots showing the ",
                ifelse(isTRUE(case$frac), "number/faction", "number"),
                " of cells per cluster",
                ifelse(
                    is.null(case$group.by),
                    "",
                    paste0(", by ", paste0(case$group.by, collapse = ", "))
                )
            )
        ),
        h1 = name
    )

    add_report(
        list(
            name = "Bar Plot",
            contents = list(list(kind = "image", src = figfile))
        ),
        h1 = name,
        ui = "tabs"
    )
    if (isTRUE(case$table)) {
        add_report(
            list(
                name = "Table",
                contents = list(list(kind = "table", src = tablefile))
            ),
            h1 = name,
            ui = "tabs"
        )
    }

    if (isTRUE(case$pie)) {
        p_pie = plot_df %>%
            arrange(!!sym(case$ident)) %>%
            ggplot(aes(x="", y=.n, fill=!!sym(case$ident))) +
            geom_bar(stat="identity", width=1, alpha=.8, position = position_stack(reverse = TRUE)) +
            coord_polar("y", start=0) +
            scale_fill_biopipen() +
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

        png(piefile, width=case$pie_devpars$width, height=case$pie_devpars$height, res=case$pie_devpars$res)
        print(p_pie)
        dev.off()

        add_report(
            list(
                name = "Pie Chart",
                contents = list(list(kind = "image", src = piefile))
            ),
            h1 = name,
            ui = "tabs"
        )
    }

    if (isTRUE(case$circos)) {
        if (is.null(case$group.by)) {
            stop(paste0(name, ": circos plots require a group-by"))
        }
        if (isTRUE(case$transpose)) {
            circos_df <- plot_df %>%
                select(from=!!sym(case$ident), to=!!sym(case$group.by), value=.n)
        } else {
            circos_df <- plot_df %>%
                select(from=!!sym(case$group.by), to=!!sym(case$ident), value=.n)
        }
        groups <- if (is.factor(circos_df$from)) {
            levels(circos_df$from)
        } else {
            unique(circos_df$from)
        }
        idents <- if (is.factor(circos_df$to)) {
            levels(circos_df$to)
        } else {
            unique(circos_df$to)
        }
        grid_cols <- pal_biopipen()(length(idents))
        names(grid_cols) <- idents
        gcols <- rep("#565656", length(groups))
        names(gcols) <- groups
        grid_cols <- c(grid_cols, gcols)
        link_cols <- grid_cols[circos_df$to]

        png(
            circosfile,
            width=case$circos_devpars$width,
            height=case$circos_devpars$height,
            res=case$circos_devpars$res
        )
        circos.clear()
        if (!isTRUE(case$circos_labels_rot)) {
            chordDiagram(
                circos_df,
                grid.col = grid_cols,
                col = link_cols,
                direction = 1,
                direction.type = c("diffHeight", "arrows"),
                link.arr.type = "big.arrow"
            )
        } else {
            chordDiagram(
                circos_df,
                grid.col = grid_cols,
                col = link_cols,
                direction = 1,
                annotationTrack = "grid",
                direction.type = c("diffHeight", "arrows"),
                link.arr.type = "big.arrow",
                preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(circos_df)))))
            )
            circos.track(track.index = 1, panel.fun = function(x, y) {
                circos.text(
                    CELL_META$xcenter, CELL_META$ylim[1],
                    CELL_META$sector.index,
                    facing = "clockwise",
                    niceFacing = TRUE,
                    adj = c(0, 0.5))
            }, bg.border = NA) # here set bg.border to NA is important
        }
        dev.off()

        add_report(
            list(
                name = "Circos plot",
                contents = list(list(kind = "image", src = circosfile))
            ),
            h1 = name,
            ui = "tabs"
        )
    }
}

sapply(names(stats), do_one_stats)
