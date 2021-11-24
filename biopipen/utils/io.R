read.table.opts = function(file, opts) {
    rncol = NULL
    if (!is.null(opts$row.names) && opts$row.names < 0) {
        rncol = -opts$row.names
        opts$row.names = NULL
        opts = c(opts, list(row.names=NULL))
    }
    if (endsWith(file, ".gz")) {
        opts$file = gzfile(file)
    } else {
        opts$file = file
    }
    out = do.call(read.table, opts)
    if (!is.null(rncol)) {
        rnames = make.unique(out[, rncol])
        out = out[, -rncol, drop=F]
        rownames(out) = rnames
    }
    return (out)
}
