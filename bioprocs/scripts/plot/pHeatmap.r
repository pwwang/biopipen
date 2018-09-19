{{rimport}}('plot.r')

library(ggplot2) # make aes available

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
cnames  = {{args.cnames | R }}
rnames  = {{args.rnames | R }}
devpars = {{args.devpars | R}}

data = read.table(infile, header = cnames, row.names = if (rnames) 1 else NULL, sep = '\t', check.names = F)

eval(parse(text = {{args.helper | repr}}))
params  = {{args.params | R}}
ggs     = {{args.ggs | R}}

plot.heatmap(data, outfile, params, ggs, devpars = devpars)
