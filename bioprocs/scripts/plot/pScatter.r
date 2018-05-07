{{rimport}}('plot.r')

infile  = {{in.infile | R}}
outfile = {{out.outfile | R}}
cnames  = {{args.cnames | R }}
rnames  = {{args.rnames | R }}
devpars = {{args.devpars | R}}
params  = {{args.params | R}}
x       = {{args.x | repr}}
y       = {{args.y | repr}}

data = read.table(infile, header = cnames, row.names = if (rnames) 1 else NULL, sep = '\t', check.names = F)
if (!is.na(as.numeric(x))) x = as.integer(x)
if (!is.na(as.numeric(y))) y = as.integer(y)

# see http://www.sthda.com/english/rpkgs/ggpubr/reference/ggscatter.html for all params
scatter(data, outfile, x, y, params, devpars = devpars)
