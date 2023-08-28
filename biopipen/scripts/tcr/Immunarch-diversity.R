# Diversity estimation
# https://immunarch.com/articles/web_only/v6_diversity.html

# Other variables are loaded in the parent template
# immdata is already loaded, meta is mutated
div_filter = {{envs.divs.filter | r}}
div_method = {{envs.divs.method | r}}
div_by = {{envs.divs.by | r}}
div_order = {{envs.divs.order | r}}
div_args = {{envs.divs.args | r: todot="-"}}
div_test = {{envs.divs.test | r}}
div_cases = {{envs.divs.cases | r: todot="-"}}
div_devpars = {{envs.divs.devpars | r}}
div_separate_by = {{envs.divs.separate_by | r}}
div_align_x = {{envs.divs.align_x | r}}
div_align_y = {{envs.divs.align_y | r}}
div_log = {{envs.divs.log | r}}

div_dir = file.path(outdir, "diversity")
dir.create(div_dir, showWarnings = FALSE)

print("- Diversity estimation ...")

# Fill up the cases
update_case = function(case) {
    if (is.null(case$filter)) {
        case$filter = div_filter
    }
    if (is.null(case$method)) {
        case$method = div_method
    }
    if (is.null(case$by)) {
        case$by = div_by
    }
    if (!is.null(case$by) && nchar(case$by) > 0) {
        case$by = unlist(strsplit(case$by, ",")) %>% trimws()
    }
    if (is.null(case$order)) {
        case$order = div_order
    }
    if (is.null(case$args)) {
        case$args = div_args
    }
    for (name in names(case$args)) {
        if (is.null(case$args[[name]])) {
            case$args[[name]] = div_args[[name]]
        }
    }
    if (is.null(case$test)) {
        case$test = div_test
    }
    if (is.null(case$test$method)) {
        case$test$method = div_test$method
    }
    if (!case$test$method %in% c("none", "t.test", "wilcox.test")) {
        stop(paste0(
            "Diversity estimation: Unknown test method: ",
            case$test$method,
            ". Expected: t.test, wilcox.test"
        ))
    }
    if (is.null(case$test$padjust)) {
        case$test$group = div_test$padjust
    }
    if (is.null(case$devpars)) {
        case$devpars = div_devpars
    }
    if (is.null(case$devpars$res)) {
        case$devpars$res = div_devpars$res
    }
    if (is.null(case$devpars$width)) {
        case$devpars$width = div_devpars$width
    }
    if (is.null(case$devpars$height)) {
        case$devpars$height = div_devpars$height
    }
    if (is.null(case$separate_by)) {
        case$separate_by = div_separate_by
    }
    if (is.null(case$align_x)) {
        case$align_x = div_align_x
    }
    if (is.null(case$align_y)) {
        case$align_y = div_align_y
    }
    if (is.null(case$log)) {
        case$log = div_log
    }
    if (!is.null(case$args) && length(case$args) > 0) {
        names(case$args) = paste0(".", names(case$args))
    }
    if (!is.null(case$test) && case$test != "none" && (is.null(case$by) || length(case$by) == 0)) {
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
    default_case = update_case(list())
    div_cases = list(x = default_case)
    names(div_cases) = div_method
} else {
    for (name in names(div_cases)) {
        div_cases[[name]] = update_case(div_cases[[name]])
    }
}

# Run different diversity estimation methods
# For general cases
# Args:
#   data: the data to be used
#   case: the case with argument to be run
#   ddir: the directory to save the results
#   trans_div: the transformation to be applied to the diversity matrix
#     note that the transformation must include the case$by columns
#   plotting: the plotting function to be used (if not "auto")
#   plot_modifier: A function to modify the plot before saving it
#   return_p: whether to return the plot
run_general = function(
    data,
    case,
    ddir,
    value_col = "Value",
    trans_div = "auto",
    plotting = "auto"
) {
    args = case$args
    args$.data = data
    args$.method = case$method
    div = do_call(repDiversity, args)
    if (is.character(trans_div) && trans_div == "auto") {
        # let's see if div has Sample column, otherwise, it should have rownames
        # as Sample
        newdiv = as.data.frame(div)
        if (!"Sample" %in% colnames(newdiv)) {
            newdiv = newdiv %>% rownames_to_column("Sample")
        }
        if (!is.null(case$by) && length(case$by) > 0) {
            newdiv = newdiv %>% left_join(
                immdata$meta[, c("Sample", case$by)],
                by = "Sample",
                suffix = c(".div", "")
            )
        }
    } else if (is.function(trans_div)) {
        newdiv = trans_div(div)
    }

    write.table(
        newdiv,
        file = file.path(ddir, "diversity.txt"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE
    )

    # plot
    #  by, order
    if (is.character(plotting) && plotting == "auto") {
        if (!is.null(case$by) && length(case$by) > 0) {
            p = vis(div, .by = case$by, .meta = immdata$meta)
        } else {
            p = vis(div)
        }
        if (!is.null(case$order) && length(case$order) > 0) {
            p = p + scale_x_discrete(limits = case$order)
        }
    } else {
        p = plotting(newdiv)
    }
    # calculate the width, height and res if not specified
    width = case$devpars$width
    height = case$devpars$height
    res = case$devpars$res
    if (is.null(res)) { res = 100 }
    if (is.null(height)) { height = 1000 }
    if (is.null(width)) {
        if (!is.null(case$by) && length(case$by) > 0) {
            width = 200 * length(unique(immdata$meta[[case$by]])) + 120
        } else {
            width = 100 * length(unique(immdata$meta$Sample)) + 120
        }
    }
    png(file.path(ddir, "diversity.png"), width = width, height = height, res = res)
    print(p)
    dev.off()

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
        tested = testfun(newdiv[[value_col]], newdiv[[groupname]], p.adjust.method = if (is.null(case$test$padjust)) "none" else case$test$padjust)
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
        write.table(
            pv,
            file = file.path(ddir, paste0("diversity.test.", case$test$method, ".txt")),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE
        )
    }
}

# gini plotting is not supported by vis()
run_gini = function(data, case, ddir) {
    trans_div = function(df) {
        df = df %>% as.data.frame() %>% rownames_to_column("Sample") %>% rename(Value = V1)

        groupname = case$by
        if (!is.null(groupname) && length(groupname) > 0) {
            df = df %>% left_join(
                immdata$meta[, c("Sample", groupname)],
                by = "Sample",
                suffix = c(".div", "")
            )
            if (length(groupname) > 1) {
                groupname = paste(groupname, collapse = "_")
                newdiv = newdiv %>% mutate(!!sym(groupname) := paste(!!!syms(case$by), collapse = "_"))
            }
        }
        return (df)
    }

    plotting = function(df) {
        if (!is.null(case$by) && length(case$by) > 0) {
            groupname = paste(case$by, collapse = "_")
            p = ggplot(df) +
                geom_boxplot(aes(x=!!sym(groupname), y=Value, fill=!!sym(groupname)), alpha = 0.5) +
                geom_jitter(aes(x=!!sym(groupname), y=Value), color="black", width=0.2, alpha=0.5) +
                theme(
                    legend.position="none",
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
                )
        } else {
            p = ggplot(df) +
                geom_col(aes(x=Sample, y=Value, fill=Sample)) +
                ggtitle("Gini Coefficient") +
                theme(
                    legend.position="none",
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
                )
        }
        if (!is.null(case$order) && length(case$order) > 0) {
            p = p + scale_x_discrete(limits = case$order)
        }
        return (p)
    }

    run_general(data, case, ddir, trans_div = trans_div, plotting = plotting)
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
        for (sam in names(idata$data)) {
            args$.data = idata$data[sam]
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
        args$.data = idata$data[valid_samples]
        do_call(repDiversity, args)
    })
}

plot_raref = function(df, log) {
    p = if (isTRUE(log)) vis(df, .log = TRUE) else vis(df)
    p + xlab("Sample size (cells)")
}

run_raref_single = function(data, case, ddir, suffix = "", save_p = TRUE) {
    args = case$args
    args$.data = data
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
        devpars$width =  750 + ceiling(length(data) / 20) * 250
    }
    if (save_p) {
        png(file.path(ddir, "raref.png"), width = devpars$width, height = devpars$height, res = devpars$res)
        print(p)
        dev.off()
    } else {
        return (list(p = p, width = devpars$width))
    }
}

run_raref_multi = function(data, case, ddir) {
    maxx = 0
    maxy = 0
    res = if (is.null(case$devpars$res)) 100 else case$devpars$res
    height = if (is.null(case$devpars$height)) 1000 else case$devpars$height
    sepvars = unique(immdata$meta[[case$separate_by]])
    width = 0
    widths = list()
    plots = list()
    for (sepvar in sepvars) {
        print(paste0("  ", case$separate_by, ": ", sepvar))
        q = list(include(sepvar))
        names(q) = case$separate_by
        single_run = run_raref_single(
            repFilter(list(data=data, meta=immdata$meta), .method = "by.meta", .query = q, .match = "exact"),
            case,
            ddir,
            suffix = paste0("-", sepvar),
            save_p = FALSE
        )
        plots[[sepvar]] = single_run$p
        if (case$align_x) {
            maxx = max(maxx, layer_scales(single_run$p)$x$range$range[2])
            width = max(width, single_run$width)
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
        } else {
            width = widths[[sepvar]]
        }
        if (case$align_y) {
            plots[[sepvar]] = plots[[sepvar]] + ylim(c(0, maxy))
        }
        png(file.path(ddir, paste0("raref-", sepvar, ".png")), width = width, height = height, res = res)
        print(plots[[sepvar]])
        dev.off()
    }
}

# Run the diversity estimation for one case
run_div_case = function(casename) {
    print(paste0("  Case: ", casename))
    case = div_cases[[casename]]
    if (case$method == "raref") {
        ddir = file.path(outdir, "rarefraction", casename)
    } else {
        ddir = file.path(div_dir, casename)
    }
    dir.create(ddir, recursive = TRUE, showWarnings = FALSE)

    # Filter
    data = immdata$data
    if (!is.null(case$filter) && length(case$filter) > 0 && nchar(case$filter) > 0) {
        for (n in names(data)) {
            data[[n]] = data[[n]] %>% filter(!!parse_expr(case$filter))
        }
    }

    # Run repDiversity
    if (case$method == "chao1") {
        run_general(data, case, ddir, "Estimator")
    } else if (case$method == "hill") {
        run_general(data, case, ddir)
    } else if (case$method == "div") {
        run_general(data, case, ddir)
    } else if (case$method == "gini.simp") {
        run_general(data, case, ddir)
    } else if (case$method == "inv.simp") {
        run_general(data, case, ddir)
    } else if (case$method == "gini") {
        run_gini(data, case, ddir)
    } else if (case$method == "raref") {
        if (!is.null(case$separate_by)) {
            run_raref_multi(data, case, ddir)
        } else {
            run_raref_single(data, case, ddir)
        }
    } else {
        stop(paste0("Unknown diversity method: ", case$method))
    }
}

# Run all cases
casenames = names(div_cases)
sapply(casenames, run_div_case)
