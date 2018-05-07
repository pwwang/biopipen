(function(...) {
	library(reticulate)
	bioprocs = import('bioprocs')
	for (rfile in list(...)) {
		source(file.path(bioprocs$UTILS, rfile))
	}
})('plot.r')

require('jsonlite')
# plot.scatter.r datafile plotfile x y params
args = commandArgs(trailingOnly=TRUE)
data = read.table(args[1], head = T, row.names = NULL, sep = '\t', check.names = F)
x    = args[3]
y    = args[4]
if (!is.na(as.numeric(x))) x = as.integer(x)
if (!is.na(as.numeric(y))) y = as.integer(y)
tryCatch({
	params = fromJSON(args[5])
}, error = function(){
	params = list()
})
scatter(data, args[2], x, y, params, devpars = list(res = 48, width = 200, height = 200))
