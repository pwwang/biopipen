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
data = read.table(args[1], head = T, row.names = 1, sep = '\t', check.names = F)

params = list()
ggs = list()
tryCatch({
	params = fromJSON(args[3])
}, error = function(e){
	params = list()
})

tryCatch({
	ggs = fromJSON(args[4])
}, error = function(e){
	ggs = list()
})
#	ggs    = list()
plot.heatmap(data, args[2], params, ggs, devpars = list(res = 48, width = 200, height = 200))
