library(ggplot2)
pdf(NULL) # preventing Rplots.pdf

plotVenn = function(
    # A named list with elements,
    # e.g. list(A=paste0("R", 1:5), B=paste0("R": 3:7))
    data,
    # Arguments for `ggVennDiagram()`
    args = list(),
    # Extra ggplot components in string
    ggs = NULL,
    # Parameters for device (res, width, height) for `png()`
    devpars = NULL,
    # The output file. If NULL, will return the plot object
    outfile = NULL
) {
    library(ggVennDiagram)

    args$x = data
    p = do.call(ggVennDiagram, args)
    if (!is.null(ggs)) {
        for (gg in ggs) {
            p = p + eval(parse(text=gg))
        }
    }

    if (is.null(outfile)) {
        return (p)
    } else {
        devpars$filename = outfile
        do.call(png, devpars)
        print(p)
        dev.off()
    }
}


plotGG = function(
    # A data frame (long format)
    data,
    # the geom
    geom,
    # Arguments for `geom_x()`
    args = list(),
    # Extra ggplot components in string
    ggs = NULL,
    # Parameters for device (res, width, height) for `png()`
    devpars = NULL,
    # The output file. If NULL, will return the plot object
    outfile = NULL
) {

    p = ggplot(data)
    p = p + do.call(paste0("geom_", geom), args)
    if (!is.null(ggs)) {
        for (gg in ggs) {
            p = p + eval(parse(text=gg))
        }
    }

    if (is.null(outfile)) {
        return (p)
    } else {
        devpars$filename = outfile
        do.call(png, devpars)
        print(p)
        dev.off()
    }
}


plotViolin = function(
    # A data frame (long format)
    data,
    # Arguments for `geom_violin()`
    args = list(),
    # Extra ggplot components in string
    ggs = NULL,
    # Parameters for device (res, width, height) for `png()`
    devpars = NULL,
    # The output file. If NULL, will return the plot object
    outfile = NULL
) {
    plotGG(data, "violin", args, ggs, devpars, outfile)
}

plotHeatmap = function(
    # Data matrix
	data,
    # Arguments for `ComplexHeatmap::Heatmap()`
    args = list(),
    # Other arguments for `ComplexHeatmap::draw()`
    draw = list(),
    # Parameters for device (res, width, height) for `png()`
    devpars = NULL,
    # The output file. If NULL, will return the plot object
    # If "draw", will call `ComplexHeatmap::draw()`
    outfile = NULL
) {
	library(ComplexHeatmap)

	args$matrix = as.matrix(data)
	hm = do.call(Heatmap, args)

    if (is.null(outfile)) {
        return(hm)
    } else if (outfile == "draw") {
        do.call(ComplexHeatmap::draw, c(list(hm), draw))
    } else {
        devpars$filename = outfile
        do.call(png, devpars)
        do.call(ComplexHeatmap::draw, c(list(hm), draw))
        dev.off()
    }
}
