{{rimport}}('__init__.r', 'plot.r')

#library(ggplot2) # make aes available

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
inopts  = {{args.inopts | R}}
devpars = {{args.devpars | R}}

data = read.table.inopts(infile, inopts)

eval(parse(text = {{args.helper | repr}}))
params  = {{args.params | R}}
ggs     = {{args.ggs | R}}

plot.heatmap(data, outfile, params, ggs, devpars = devpars)
