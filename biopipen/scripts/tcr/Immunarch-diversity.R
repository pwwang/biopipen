
# https://immunarch.com/articles/web_only/v6_diversity.html

log_info("")
log_info("# Diversity estimation")
log_info("-----------------------------------")

div_test = div_test %||% list(method = "none", padjust = "none")
div_devpars = div_devpars %||% list(res = 100, width = 800, height = 800)

div_dir = file.path(outdir, "diversity")
dir.create(div_dir, showWarnings = FALSE)

# Fill up the cases
update_case = function(case, name) {
    log_debug("Filling up case: {name} ...")
    case$subset <- case$subset %||% div_subset
    case$method <- case$method %||% div_method
    case$by <- case$by %||% div_by
    if (!is.null(case$by) && nchar(case$by) > 0) {
        case$by = unlist(strsplit(case$by, ",")) %>% trimws()
    }
    case$plot_type <- case$plot_type %||% div_plot_type
    case$order <- case$order %||% div_order
    case$args <- case$args %||% div_args
    for (name in names(case$args)) {
        case$args[[name]] = case$args[[name]] %||% div_args[[name]]
    }
    case$test <- case$test %||% div_test
    case$test$method <- case$test$method %||% div_test$method
    if (!case$test$method %in% c("none", "t.test", "wilcox.test")) {
        stop(paste0(
            "Diversity estimation: Unknown test method: ",
            case$test$method,
            ". Expected: t.test, wilcox.test"
        ))
    }
    case$test$padjust <- case$test$padjust %||% div_test$padjust
    case$devpars <- case$devpars %||% div_devpars
    case$devpars$res <- case$devpars$res %||% div_devpars$res
    case$devpars$width <- case$devpars$width %||% div_devpars$width
    case$devpars$height <- case$devpars$height %||% div_devpars$height
    case$separate_by <- case$separate_by %||% div_separate_by
    case$split_by <- case$split_by %||% div_split_by
    case$split_order <- case$split_order %||% div_split_order
    if (!is.null(case$separate_by) && !is.null(case$split_by)) {
        stop("Diversity estimation: `separate_by` and `split_by` cannot be specified at the same time")
    }
    case$align_x <- case$align_x %||% div_align_x
    case$align_y <- case$align_y %||% div_align_y
    case$log <- case$log %||% div_log
    case$ncol <- case$ncol %||% div_ncol
    case$ymin <- case$ymin %||% div_ymin
    case$ymax <- case$ymax %||% div_ymax
    if (!is.null(case$args) && length(case$args) > 0) {
        names(case$args) = paste0(".", names(case$args))
    }
    if (!is.null(case$test) && case$test$method != "none" && (is.null(case$by) || length(case$by) == 0)) {
        stop("For diversity estimation, `test` is only supported when `by` is specified")
    }
    # Just ignore them for rarefraction
    # if (!is.null(case$test) && case$method == "raref") {
    #     stop("For diversity estimation, `test` is not supported when `method` is `raref`")
    # }
    # if (!is.null(case$order) && case$method == "raref") {
    #     stop("For diversity estimation, `order` is not supported when `method` is `raref`")
    # }
    return (case)
}

if (is.null(div_cases) || length(div_cases) == 0) {
    if (is.null(div_method) || length(div_method) == 0 || nchar(div_method) == 0) {
        stop("No method is specified for diversity estimation")
    }
    default_case = update_case(list(), name = "DEFAULT")
    div_cases = list(x = default_case)
    names(div_cases) = div_method
} else {
    for (name in names(div_cases)) {
        div_cases[[name]] = update_case(div_cases[[name]], name = name)
    }
}

filter_div = function(div, samples) {
    if ("Sample" %in% colnames(div)) {
        dv = div[div$Sample %in% samples, , drop = FALSE]
    } else {
        dv = div[rownames(div) %in% samples, , drop = FALSE]
    }
    class(dv) = class(div)
    dv
}

# Run different diversity estimation methods
# For general cases
# Args:
#   d: the data to be used
#   case: the case with argument to be run
#   ddir: the directory to save the results
#   value_col: the column name of the value
run_general = function(casename, d, case, ddir, value_col = "Value") {
    args = case$args
    args$.data = d$data
    args$.method = case$method
    div = do_call(repDiversity, args)
    # let's see if div has Sample column, otherwise, it should have rownames
    # as Sample
    newdiv = as.data.frame(div)
    if (!"Sample" %in% colnames(newdiv)) {
        newdiv = newdiv %>% rownames_to_column("Sample")
    }
    if (!is.null(case$by) && length(case$by) > 0) {
        newdiv = newdiv %>% left_join(
            d$meta[, c("Sample", case$by), drop = FALSE],
            by = "Sample",
            suffix = c(".div", "")
        )
    }

    if (!is.null(case$separate_by)) {
        newdiv = newdiv %>% left_join(
            d$meta[, c("Sample", case$separate_by), drop = FALSE],
            by = "Sample",
            suffix = c(".div", "")
        )
    }

    if (!is.null(case$split_by)) {
        newdiv = newdiv %>% left_join(
            d$meta[, c("Sample", case$split_by), drop = FALSE],
            by = "Sample",
            suffix = c(".div", "")
        )
    }

    write.table(
        newdiv,
        file = file.path(ddir, "diversity.txt"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE
    )

    .meta_vals <- function(meta, cols) {
        if (length(cols) == 1) {
            return (meta[[cols]])
        }

        vlist = lapply(cols, function(.x) meta[[.x]])
        do.call(function(...) paste(..., sep = "; "), vlist)
    }

    # plot
    #  by, order, separate_by, align_y
    n_seps = 1
    if (!is.null(case$by) && length(case$by) > 0) {
        if (!is.null(case$separate_by)) {
            metas = split(d$meta, d$meta[[case$separate_by]])
            if (!is.null(case$split_order)) {
                if (is.character(case$split_order) && length(case$split_order) == 1) {
                    case$split_order = trimws(unlist(strsplit(case$split_order, ",")))
                }
                metas = metas[intersect(case$split_order, names(metas))]
            }
            ps = lapply(metas, function(meta) {
                .test = length(unique(.meta_vals(meta, case$by))) > 1
                p = vis(
                    filter_div(div, meta$Sample),
                    .by = case$by,
                    .meta = meta,
                    .test = .test,
                    .plot.type = case$plot_type
                )
                p = p + xlab(paste0(case$separate_by, ": ", meta[[case$separate_by]][1], ")"))
                if (!is.null(case$order) && length(case$order) > 0) {
                    p = p + scale_x_discrete(
                        limits = intersect(case$order, unique(.meta_vals(meta, case$by)))
                    )
                }
                if (!is.null(case$ymin) && !is.null(case$ymax)) {
                    p = p + ylim(c(case$ymin, case$ymax))
                } else if (case$align_y) {
                    m1 = min(newdiv[[value_col]])
                    m2 = max(newdiv[[value_col]])
                    margin = (m2 - m1) * 0.1
                    p = p + ylim(c(m1 - margin, m2 + margin))
                }
                p
            })
            n_seps = length(ps)
            p = wrap_plots(ps, ncol = case$ncol)
        } else if (!is.null(case$split_by)) {
            metas = split(d$meta, d$meta[[case$split_by]])
            if (!is.null(case$split_order)) {
                if (is.character(case$split_order) && length(case$split_order) == 1) {
                    case$split_order = trimws(unlist(strsplit(case$split_order, ",")))
                }
                metas = metas[intersect(case$split_order, names(metas))]
            }
            .i = 0
            ps = lapply(metas, function(meta) {
                nby = length(unique(.meta_vals(meta, case$by)))
                p = vis(
                    filter_div(div, meta$Sample),
                    .by = case$by,
                    .meta = meta,
                    .test = nby > 1,
                    .plot.type = case$plot_type
                )
                if (!is.null(case$order) && length(case$order) > 0) {
                    p = p + scale_x_discrete(
                        limits = intersect(case$order, unique(.meta_vals(meta, case$by)))
                    )
                }
                p = p + xlab(meta[[case$split_by]][1]) + theme(
                    axis.text.x = element_blank(),
                    plot.title = element_blank(),
                    plot.subtitle = element_blank(),
                    legend.position = "right"
                )
                if (.i > 0) {
                    p = p + theme(
                        axis.text.y = element_blank(),
                        axis.title.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.line.y = element_blank()
                    )
                }
                .i <<- .i + 1
                list(
                    p = p,
                    ymin = case$ymin %||% min(newdiv[[value_col]]),
                    ymax = case$ymax %||% max(newdiv[[value_col]]),
                    nby = nby
                )
            })
            n_seps = length(ps)
            ymin = do.call(min, lapply(ps, function(x) x$ymin))
            ymin = max(0, ymin - 0.1)
            ymax = do.call(max, lapply(ps, function(x) x$ymax))
            ymax = ymax + 0.3 * (ymax - ymin)  # for the pvalue marks
            widths = sapply(ps, function(x) ifelse(x$nby == 1, 1.2, x$nby))
            plots = lapply(ps, function(x) x$p + ylim(c(ymin, ymax)))
            p = wrap_plots(plots, widths = widths, guides = "collect")
        } else {
            .test = length(unique(.meta_vals(d$meta, case$by))) > 1
            p = vis(div, .by = case$by, .meta = d$meta, .test = .test, .plot.type = case$plot_type)
            if (!is.null(case$order) && length(case$order) > 0) {
                p = p + scale_x_discrete(limits = intersect(case$order, unique(.meta_vals(d$meta, case$by))))
            }
        }
    } else if (!is.null(case$separate_by)) {
        metas = split(d$meta, d$meta[[case$separate_by]])
        if (!is.null(case$split_order)) {
            if (is.character(case$split_order) && length(case$split_order) == 1) {
                case$split_order = trimws(unlist(strsplit(case$split_order, ",")))
            }
            metas = metas[intersect(case$split_order, names(metas))]
        }
        ps = lapply(metas, function(meta) {
            p = vis(filter_div(div, meta$Sample))
            p = p + ggtitle(paste0(p$labels$title, " (" , case$separate_by, ": ", meta[[case$separate_by]][1], ")"))
            if (!is.null(case$order) && length(case$order) > 0) {
                p = p + scale_x_discrete(limits = intersect(case$order, unique(meta[[Sample]])))
            }
            if (!is.null(case$ymin) && !is.null(case$ymax)) {
                p = p + ylim(c(case$ymin, case$ymax))
            } else if (case$align_y) {
                m1 = min(newdiv[[value_col]])
                m2 = max(newdiv[[value_col]])
                margin = (m2 - m1) * 0.1
                p = p + ylim(c(m1 - margin, m2 + margin))
            }
            p
        })
        n_seps = length(ps)
        p = wrap_plots(ps, ncol = case$ncol)
    } else if (!is.null(case$split_by)) {
        metas = split(d$meta, d$meta[[case$split_by]])
        if (!is.null(case$split_order)) {
            if (is.character(case$split_order) && length(case$split_order) == 1) {
                case$split_order = trimws(unlist(strsplit(case$split_order, ",")))
            }
            metas = metas[intersect(case$split_order, names(metas))]
        }
        .i = 0
        ps = lapply(metas, function(meta) {
            nby = length(unique(meta$Sample))
            p = vis(filter_div(div, meta$Sample))
            if (!is.null(case$order) && length(case$order) > 0) {
                p = p + scale_x_discrete(limits = intersect(case$order, unique(meta[[Sample]])))
            }
            p = p + xlab(meta[[case$split_by]][1]) + theme(
                axis.text.x = element_blank(),
                plot.title = element_blank(),
                plot.subtitle = element_blank(),
                legend.position = "right"
            )
            if (.i > 0) {
                p = p + theme(
                    axis.text.y = element_blank(),
                    axis.title.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.line.y = element_blank()
                )
            }
            .i <<- .i + 1
            list(
                p = p,
                ymin = case$ymin %||% min(newdiv[[value_col]]),
                ymax = case$ymax %||% max(newdiv[[value_col]]) + 0.1 * max(newdiv[[value_col]]),
                nby = nby
            )
        })
        n_seps = length(ps)
        ymin = do.call(min, lapply(ps, function(x) x$ymin))
        ymin = max(0, ymin - 0.1)
        ymax = do.call(max, lapply(ps, function(x) x$ymax))
        ymax = ymax + 0.3 * (ymax - ymin)  # for the pvalue marks
        widths = sapply(ps, function(x) ifelse(x$nby == 1, 1.2, x$nby))
        plots = lapply(ps, function(x) x$p + ylim(c(ymin, ymax)))
        p = wrap_plots(plots, widths = widths, guides = "collect")
    } else {
        p = vis(div)
        if (!is.null(case$order) && length(case$order) > 0) {
            p = p + scale_x_discrete(
                limits = intersect(case$order, unique(.meta_vals(d$meta, case$by)))
            )
        }
    }

    # calculate the width, height and res if not specified
    width = case$devpars$width
    height = case$devpars$height
    res = case$devpars$res
    res = res %||% 100
    if (is.null(height)) {
        if (!is.nulL(case$split_by)) {
            height = 800
        } else {
            height = if (n_seps == 1) 800 else 600 * ceiling(n_seps / case$ncol)
        }
    }
    if (is.null(width)) {
        if (!is.null(case$by) && length(case$by) > 0) {
            width = 200 * length(unique(.meta_vals(d$meta, case$by))) + 120
        } else {
            width = 100 * length(unique(d$meta$Sample)) + 120
        }
        if (!is.null(case$split_by)) { width = width / 2 }
        if (n_seps > 1) {
            width = width * case$ncol
        }
    }

    div_plot = file.path(ddir, "diversity.png")
    png(div_plot, width = width, height = height, res = res)
    print(p)
    dev.off()

    div_plot_pdf = file.path(ddir, "diversity.pdf")
    pdf(div_plot_pdf, width = width / res, height = height / res)
    print(p)
    dev.off()

    add_report(
        list(
            kind = "descr",
            content = paste0(
                "Diversity estimation using ",
                "<code>",
                case$method,
                "</code>: ",
                switch(case$method,
                    chao1 = paste0(
                        "a nonparameteric asymptotic estimator of species richness ",
                        "(number of species in a population)."),
                    hill = paste0(
                        "Hill numbers are a mathematically unified family of ",
                        "diversity indices (differing only by an exponent q)."),
                    div = paste0(
                        "true diversity, or the effective number of types, ",
                        "refers to the number of equally abundant types needed for ",
                        "the average proportional abundance of the types to equal that ",
                        "observed in the dataset of interest where all types may ",
                        "not be equally abundant."),
                    gini.simp = paste0(
                        "the Gini-Simpson index is the probability of interspecific ",
                        "encounter, i.e., probability that two entities represent different types."),
                    inv.simp = paste0(
                        "Inverse Simpson index is the effective number of types ",
                        "that is obtained when the weighted arithmetic mean is used ",
                        "to quantify average proportional abundance of types in ",
                        "the dataset of interest."),
                    gini = paste0(
                        "the Gini coefficient measures the inequality among ",
                        "values of a frequency distribution (for example levels of income). ",
                        "A Gini coefficient of zero expresses perfect equality, ",
                        "where all values are the same (for example, where everyone has ",
                        "the same income). A Gini coefficient of one (or 100 percents ) ",
                        "expresses maximal inequality among values (for example where only ",
                        "one person has all the income)."),
                    d50 = paste0(
                        "the D50 index. ",
                        "It is the number of types that are needed to cover 50% of the total
                        abundance.")
                )
            )
        ),
        h1 = "Diversity Estimation",
        h2 = casename
    )
    add_report(
        list(
            name = "Diversity Plot",
            contents = list(list(kind = "image", src = div_plot, download = div_plot_pdf))
        ),
        list(
            name = "Diversity Table",
            contents = list(list(kind = "table", src = file.path(ddir, "diversity.txt")))
        ),
        h1 = "Diversity Estimation",
        h2 = casename,
        ui = "tabs"
    )

    # Test
    if (!is.null(case$test) && case$test$method != "none") {
        # Use pairwise.t.test or pairwise.wilcox.test
        if (case$test$method == "t.test") {
            testfun = pairwise.t.test
        } else if (case$test$method == "wilcox.test") {
            testfun = pairwise.wilcox.test
        }

        groupname = case$by
        if (length(groupname) > 1) {
            groupname = paste(groupname, collapse = "_")
            newdiv = newdiv %>% mutate(!!sym(groupname) := paste(!!!syms(case$by), collapse = "_"))
        }

        get_pv = function(nd) {
            if (length(unique(nd[[groupname]])) <= 1) {
                return (NULL)
            }
            tested = testfun(nd[[value_col]], nd[[groupname]], p.adjust.method = if (is.null(case$test$padjust)) "none" else case$test$padjust)
            pv = as.data.frame(tested$p.value)
            # Make pv symmetric
            #  From
            #         A            B
            # B 0.4322773         NA
            # C 0.1088602 0.08353865
            #  To
            #         A            B            C
            # A        NA 0.432277300 0.108860200
            # B 0.4322773           NA 0.083538650
            # C 0.1088602 0.083538650           NA
            rname_to_add = colnames(pv)[1]
            cname_to_add = rownames(pv)[nrow(pv)]
            pv[[cname_to_add]] = pv[nrow(pv), ] %>% unlist() %>% unname()
            pv = rbind(NA, pv)
            rownames(pv)[1] = rname_to_add
            pv[[cname_to_add]] = shift(pv[[cname_to_add]], n = -1, fill = NA)
            pv[upper.tri(pv)] = pv[lower.tri(pv)]
            pv = formatC(as.matrix(pv), format = "e", digits = 2) %>%
                trimws() %>%
                as.data.frame() %>%
                rownames_to_column("Group")
        }
        if (!is.null(case$separate_by)) {
            pvs = lapply(split(newdiv, newdiv[[case$separate_by]]), function(nd) {
                pv = get_pv(nd)
                pv[[case$separate_by]] = nd[[case$separate_by]][1]
                pv
            })
            pv = do.call(rbind, pvs)
        } else {
            pv = get_pv(newdiv)
        }

        write.table(
            pv,
            file = file.path(ddir, paste0("diversity.test.", case$test$method, ".txt")),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE
        )

        add_report(
            list(
                name = paste0("Test (", case$test$method, ")"),
                contents = list(list(
                    kind = "table",
                    src = file.path(ddir, paste0("diversity.test.", case$test$method, ".txt"))
                ))
            ),
            h1 = "Diversity Estimation",
            h2 = casename,
            ui = "tabs"
        )
    }
}

# rarefraction
raref_rep_diversity = function(args) {
    # Due to https://github.com/immunomind/immunarch/issues/44
    idata = args$.data

    tryCatch({
        do_call(repDiversity, args)
    }, error = function(e) {
        # https://github.com/immunomind/immunarch/issues/44
        valid_samples = c()
        for (sam in names(idata)) {
            args$.data = idata[sam]
            vsam = tryCatch({
                do_call(repDiversity, args)
                sam
            }, error=function(e) {
                warning(
                    paste(
                        "Rarefraction analysis failed for sample",
                        sam,
                        ":",
                        as.character(e)
                    ),
                    immediate. = TRUE
                )
                c()
            })
            valid_samples = c(valid_samples, vsam)
        }
        args$.data = idata[valid_samples]
        do_call(repDiversity, args)
    })
}

plot_raref = function(df, log) {
    p = if (isTRUE(log)) vis(df, .log = TRUE) else vis(df)
    p + xlab("Sample size (cells)")
}

run_raref_single = function(d, case, ddir, suffix = "", save_p = TRUE) {
    args = case$args
    args$.data = d$data
    args$.method = case$method
    div = raref_rep_diversity(args)
    write.table(
        div,
        file = file.path(ddir, paste0("raref", suffix, ".txt")),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE
    )
    p = plot_raref(div, case$log)
    devpars = case$devpars
    if (is.null(devpars$res)) { devpars$res = 100 }
    if (is.null(devpars$height)) { devpars$height = 1000 }
    if (is.null(devpars$width)) {
        devpars$width =  750 + ceiling(length(args$.data) / 20) * 250
    }
    if (save_p) {
        png(file.path(ddir, "raref.png"), width = devpars$width, height = devpars$height, res = devpars$res)
        print(p)
        dev.off()

        pdf(file.path(ddir, "raref.pdf"), width = devpars$width / devpars$res, height = devpars$height / devpars$res)
        print(p)
        dev.off()
    } else {
        return (list(p = p, width = devpars$width))
    }
}

run_raref_multi = function(d, case, ddir) {
    maxx = 0
    maxy = 0
    res = if (is.null(case$devpars$res)) 100 else case$devpars$res
    sepvars = unique(d$meta[[case$separate_by]])
    ncol = if (is.null(case$ncol)) 2 else case$ncol
    running_width = 0
    widths = list()
    plots = list()
    for (sepvar in sepvars) {
        # print(paste0("  ", case$separate_by, ": ", sepvar))
        log_info("  {case$separate_by}: {sepvar}")
        q = list(include(sepvar))
        names(q) = case$separate_by
        single_run = run_raref_single(
            repFilter(list(data=d$data, meta=d$meta), .method = "by.meta", .query = q, .match = "exact"),
            case,
            ddir,
            suffix = paste0("-", sepvar),
            save_p = FALSE
        )
        plots[[sepvar]] = single_run$p
        if (case$align_x) {
            maxx = max(maxx, layer_scales(single_run$p)$x$range$range[2])
            running_width = max(running_width, single_run$width)
        } else {
            widths[[sepvar]] = single_run$width
        }
        if (case$align_y) {
            maxy = max(maxy, layer_scales(single_run$p)$y$range$range[2])
        }
    }
    for (sepvar in sepvars) {
        if (case$align_x) {
            if (case$log) {
                plots[[sepvar]] = plots[[sepvar]] + scale_x_log10(limits = c(1, 10 ^ maxx))
            } else {
                plots[[sepvar]] = plots[[sepvar]] + xlim(c(0, maxx))
            }
        }
        if (case$align_y) {
            plots[[sepvar]] = plots[[sepvar]] + ylim(c(0, maxy))
        }
    }
    p = wrap_plots(plots, ncol = ncol)
    if (is.null(case$devpars$width)) {
        width = if (case$align_x) running_width * ncol else max(unlist(widths)) * ncol
    } else {
        width = case$devpars$width
    }
    if (is.null(case$devpars$height)) {
        nrow = ceiling(length(plots) / ncol)
        height = if (case$align_y) maxy * nrow else 600 * nrow
    } else {
        height = case$devpars$height
    }
    png(
        file.path(ddir, paste0("raref-", slugify(case$separate_by), ".png")),
        width = width,
        height = height,
        res = res
    )
    print(p)
    dev.off()

    pdf(
        file.path(ddir, paste0("raref-", slugify(case$separate_by), ".pdf")),
        width = width / res,
        height = height / res
    )
    print(p)
    dev.off()
}

# Run the diversity estimation for one case
run_div_case = function(casename) {
    log_info("Processing case: {casename} ...")
    case = div_cases[[casename]]
    if (case$method == "raref") {
        ddir = file.path(outdir, "rarefraction", slugify(casename))
    } else {
        ddir = file.path(div_dir, slugify(casename))
    }
    dir.create(ddir, recursive = TRUE, showWarnings = FALSE)

    # Filter
    if (!is.null(case$subset)) {
        d = immdata_from_expanded(filter_expanded_immdata(exdata, case$subset))
        if (nrow(d$meta) == 0) {
            stop(paste0(
                "No samples/cells left after filtering. ",
                "Do you have the correct `subset` for case: ",
                casename, "?"))
        }
    } else {
        d = immdata
    }

    # Run repDiversity
    if (case$method == "raref") {
        add_report(
            list(
                kind = "descr",
                content = paste0(
                    "Rarefaction is a technique to assess species richness from the ",
                    "results of sampling through extrapolation. "
                )
            ),
            h1 = "Rarefraction",
            h2 = casename
        )

        if (!is.null(case$separate_by)) {
            run_raref_multi(d, case, ddir)
            add_report(
                list(
                    kind = "image",
                    src = file.path(ddir, paste0("raref-", slugify(case$separate_by), ".png")),
                    download = file.path(ddir, paste0("raref-", slugify(case$separate_by), ".pdf"))
                ),
                h1 = "Rarefraction",
                h2 = casename
            )
        } else {
            run_raref_single(d, case, ddir)
            add_report(
                list(
                    kind = "image",
                    src = file.path(ddir, "raref.png"),
                    download = file.path(ddir, "raref.pdf")
                ),
                h1 = "Rarefraction",
                h2 = casename
            )
        }
    } else {
        if (case$method == "chao1") {
            run_general(casename, d, case, ddir, "Estimator")
        } else if (case$method == "hill") {
            run_general(casename, d, case, ddir)
        } else if (case$method == "div") {
            run_general(casename, d, case, ddir)
        } else if (case$method == "gini.simp") {
            run_general(casename, d, case, ddir)
        } else if (case$method == "inv.simp") {
            run_general(casename, d, case, ddir)
        } else if (case$method == "gini") {
            run_general(casename, d, case, ddir, "V1")
        } else if (case$method == "d50") {
            run_general(casename, d, case, ddir, "Clones")
        } else {
            stop(paste0("Unknown diversity method: ", case$method))
        }
    }
}

# Run all diversity cases
sapply(names(div_cases), run_div_case)
