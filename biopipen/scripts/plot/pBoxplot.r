{{rimport}}('__init__.r', 'plot.r')

library(ggplot2) # make aes available

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
inopts  = {{args.inopts | R}}
devpars = {{args.devpars | R}}
x       = {{args.x | repr}}
y       = {{args.y | repr}}
stacked = {{args.stacked | R}}

inopts$fill = TRUE
data   = read.table.inopts(infile, inopts)

if (!is.na(as.numeric(x))) x = as.integer(x)
if (!is.na(as.numeric(y))) y = as.integer(y)

eval(parse(text = {{args.helper | repr}}))
params = {{args.params | R}}
ggs    = {{args.ggs | R}}

plot.boxplot(data, outfile, x, y, stacked, params, ggs, devpars = devpars)
