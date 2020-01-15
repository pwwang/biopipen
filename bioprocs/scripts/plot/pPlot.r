{{rimport}}('plot.r')

library(ggplot2) # make aes available

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
cnames  = {{args.cnames | R }}
rnames  = {{args.rnames | R }}
devpars = {{args.devpars | R}}
aesList = {{args.aes | R}}

data = read.table(infile, header = cnames, row.names = if (rnames) 1 else NULL, sep = '\t', check.names = F)

eval(parse(text = {{args.helper | repr}}))
params  = {{args.params | R}}
ggs     = {{args.ggs | R}}

aesNames = names(aesList)
for (name in aesNames) {
	if (name == 'x' && !is.na(as.numeric(aesList$x))) {
		aesList$x = as.numeric(aesList$x)
	}
	if (name == 'y' && !is.na(as.numeric(aesList$y))) {
		aesList$y = as.numeric(aesList$y)
	}
}
if ('x' %in% aesNames && 'y' %in% aesNames) {
	plot.xy(data, outfile, x, y, params, ggs, devpars = devpars)
} else if ('x' %in% aesNames) {
	plot.x(data, outfile, x, params, ggs, devpars = devpars)
} else {
	plot.no(data, outfile, params, ggs, devpars = devpars)
}
