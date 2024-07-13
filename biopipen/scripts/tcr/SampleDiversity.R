{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
library(immunarch)

immdatafile = {{in.immdata | r}}
outputdir = {{out.outdir | r}}
# list(
#    div = c("Status", "Sex", "Status,Sex"),
#    ...
# )
div_methods = {{envs.div_methods | r}}
devpars = {{envs.devpars | r}}

immdata = readRDS(immdatafile)

# Available methods
DIV_METHODS = c("chao1", "hill", "div", "gini.simp", "inv.simp", "gini", "raref")

methods = names(div_methods)


# immunarch missing vis for immunr_gini
vis.immunr_gini <- function(.data, .by = NA, .meta = NA,
                            .errorbars = c(0.025, 0.975), .errorbars.off = FALSE,
                            .points = TRUE, .test = TRUE, .signif.label.size = 3.5, ...) {
    #         [,1]
    # s0301 0.4027197
    # s0302 0.3504762
    # s0303 0.4217427
    # s0304 0.3857278
    # s0305 0.4527641
    # attr(,"class")
    # [1] "immunr_gini" "matrix"      "array"
    .data = data.frame(Sample = rownames(.data), Value = .data[, 1])

    vis_bar(
        .data = .data, .by = .by, .meta = .meta,
        .errorbars = .errorbars, .errorbars.off = .errorbars.off, .stack = FALSE,
        .points = .points, .test = .test, .signif.label.size = .signif.label.size,
        .defgroupby = "Sample", .grouping.var = "Group",
        .labs = c(NA, "Gini Coefficient"),
        .title = "Gini Coefficient", .subtitle = "Sample diversity estimation using the Gini coefficient",
        .legend = NA, .leg.title = NA
    )
}

do_method = function(method) {
    if (!method %in% DIV_METHODS) {
        stop("Unknown method: ", method)
    } else {
        print(paste0("Method: ", method))
    }
    outdir = file.path(outputdir, method)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    # For each sample without grouping
    print("  No grouping")
    div = repDiversity(immdata$data, .method = method)

    p = vis(div)
    args = devpars[[method]]
    if (is.null(args)) {
        args = list(width = 1000, height = 1000, res = 100)
    }
    args$filename = file.path(outdir, "_nogrouping.png")
    do_call(png, args)
    print(p)
    dev.off()

    # Groupings
    bys = div_methods[[method]]
    for (by in bys) {
        print(paste0("  Grouping: ", by))
        p = tryCatch({
            vis(div, .by = unlist(strsplit(by, ",", fixed=TRUE)), .meta = immdata$meta)
        }, error = function(e) {
            print(paste("    Skipping:", e))
            return(NULL)
        })
        if (is.null(p)) {
            next
        }
        args$filename = file.path(outdir, paste0("By ", by, ".png"))
        do_call(png, args)
        print(p)
        dev.off()
    }
}

sapply(methods, do_method)
