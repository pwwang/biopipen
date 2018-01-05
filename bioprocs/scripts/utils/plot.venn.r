if (!exists('plotVenn')) {
	plotVenn = function(mat, filename, params, devpars) {
		library(VennDiagram)
		rnames = rownames(mat)
		if (is.null(rnames)) {
			rnames = paste('ROW', 1:nrow(mat), sep = '')
			rownames(mat) = rnames
		}
		x = list()
		for (cname in colnames(mat)) {
			x[[cname]] = rownames(mat[which(mat[, cname] == 1), cname, drop=F])
		}
		ps = list(x = x, filename = filename, height = devpars$height, width = devpars$width, resolution = devpars$res, imagetype = "png", alpha = .5, fill = rainbow(ncol(mat)))

		do.call(venn.diagram, c(ps, params))
		
	}
}