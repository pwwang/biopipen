{{rimport}}('plot.r')

ggs     = {{args.ggs | R}}
infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
rnames  = {{args.rnames | R}}
ggs     = {{args.ggs | R}}
devpars = {{args.devpars | R}}
data    = read.table.nodup(infile, sep = "\t", header = T, row.names = rnames)

plot.pie(data, outfile, ggs, devpars)
