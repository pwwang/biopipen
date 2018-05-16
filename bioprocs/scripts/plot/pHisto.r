{{rimport}}('plot.r')

library(ggplot2) # make aes available

infile  = {{in.infile | R}}
outfile = {{out.outfile | R}}
cnames  = {{args.cnames | R }}
rnames  = {{args.rnames | R }}
devpars = {{args.devpars | R}}
x       = {{args.x | repr}}

data = read.table(infile, header = cnames, row.names = if (rnames) 1 else NULL, sep = '\t', check.names = F)
if (!is.na(as.numeric(x))) x = as.integer(x)

eval(parse(text = {{args.helper | repr}}))
params  = {{args.params | R}}
ggs     = {{args.ggs | R}}

plot.histo(data, outfile, x, params, ggs, devpars = devpars)
