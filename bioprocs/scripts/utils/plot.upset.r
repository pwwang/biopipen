
if (!exists('plotUpset')) {
	library(UpSetR)
	plotUpset = function(mat, filename, params, devpars) {
		library(UpSetR)
		do.call(png, c(list(filename = filename), devpars))
		png (filename, res=300, width=2000, height=2000)
		if (! "nintersects" %in% names(params)) {
			params$nintersects = NA
		}
		v = do.call (upset, c(list(data = mat, sets = colnames(mat)), params))
		print (v)
		dev.off()
	}
}
