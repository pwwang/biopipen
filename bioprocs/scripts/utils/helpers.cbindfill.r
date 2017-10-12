if (!exists('cbindfill')) {
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
}