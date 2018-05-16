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

params = list()
ggs = list()
tryCatch({
	params = fromJSON(args[5])
}, error = function(e){
	params = list()
})

tryCatch({
	ggs = fromJSON(args[6])
}, error = function(e){
	ggs = list()
})
#	ggs    = list()
plot.points(data, args[2], x, y, params, ggs, devpars = list(res = 48, width = 200, height = 200))
