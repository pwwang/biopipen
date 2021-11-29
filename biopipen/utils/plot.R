library(ggplot2)

plotVenn = function(data, args, ggs, devpars = NULL, outfile = NULL) {
    # data should be a named list with elements
    library(ggVennDiagram)

    args$x = data
    p = do.call(ggVennDiagram, args)
    if (!is.null(ggs)) {
        p = p + eval(parse(text=ggs))
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