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
if (!is.na(as.numeric(x))) x = as.integer(x)

params = list()
ggs = list()
tryCatch({
	params = fromJSON(args[4])
}, error = function(e){
	params = list()
})

tryCatch({
	ggs = fromJSON(args[5])
}, error = function(e){
	ggs = list()
})
#	ggs    = list()
plot.freqpoly(data, args[2], x, params, ggs, devpars = list(res = 48, width = 200, height = 200))
