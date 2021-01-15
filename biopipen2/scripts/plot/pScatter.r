{{rimport}}('plot.r')

library(ggplot2) # make aes available

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
cnames  = {{args.cnames | R }}
rnames  = {{args.rnames | R }}
devpars = {{args.devpars | R}}
x       = {{args.x | repr}}
y       = {{args.y | repr}}

data = read.table(infile, header = cnames, row.names = if (rnames) 1 else NULL, sep = '\t', check.names = F)
if (!is.na(as.numeric(x))) x = as.integer(x)
if (!is.na(as.numeric(y))) y = as.integer(y)

eval(parse(text = {{args.helper | repr}}))
params  = {{args.params | R}}
ggs     = {{args.ggs | R}}

plot.scatter(data, outfile, x, y, params, ggs, devpars = devpars)
