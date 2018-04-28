library(reticulate)
utils     = import('bioprocs')$utils
runcmd    = utils$runcmd
mem2      = utils$mem2
# key orders not kept!
cmdargs   = utils$cmdargs

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

logger = function(..., level = 'INFO') {
	cat(paste0(level, ': ', paste(...), '\n'), file = stderr())		
}

# avoid typos
cbind.fill = cbindfill
rbind.fill = rbindfill

read.table.nodup = function(...) {
	args = list(...)
	if (!'row.names' %in% names(args) || is.null(args$row.names)) {
		return(read.table(...))
	} else {
		args$row.names = NULL
		mat = do.call(read.table, args)
		rnames = make.unique(as.character(as.vector(mat[,1,drop = T])))
		mat = mat[,-1,drop=F]
		row.names(mat)  = rnames
		return(mat)
	}
}
