if (!exists('rbindfill')) {
	rbindfill = function (x1, x2) {
		if (is.null(x1)) return(x2)
		if (is.null(x2)) return(x1)
		y = merge(x1, x2, all=T, sort=F)
		rownames(y) = c(rownames(x1), rownames(x2))
		return (y)
	}
}