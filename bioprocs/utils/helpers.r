library(reticulate)
helpers_  = import('bioprocs')$utils$helpers
runcmd    = helpers_$runcmd
mem2      = helpers_$mem2
# key orders not kept!
cmdargs   = helpers_$cmdargs

cbindfill = function (x1, x2) {
	y = merge(x1, x2, by='row.names', all=T, sort=F)
	rownames(y) = y[, "Row.names"]
	y = y[, -1, drop=F]
	cnames      = c(colnames(x1), colnames(x2))
	if (!is.null(cnames)) {
		colnames(y) = cnames
	}
	return (y)
}

rbindfill = function (x1, x2) {
	if (is.null(x1)) return(x2)
	if (is.null(x2)) return(x1)
	y = merge(x1, x2, all=T, sort=F)
	rownames(y) = make.unique(c(rownames(x1), rownames(x2)))
	return (y)
}

