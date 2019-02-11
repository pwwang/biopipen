{{rimport}}('__init__.r', 'plot.r')

library(ggplot2) # make aes available

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
inopts  = {{args.inopts | R}}
devpars = {{args.devpars | R}}
params  = {{args.params | R}}
ggs     = {{args.ggs | R}}
x       = {{args.x | R}}
y       = {{args.y | R}}
stacked = {{args.stacked | R}}

inopts$file        = infile
inopts$header      = inopts$cnames
inopts$cnames      = NULL
inopts$sep         = inopts$delimit
inopts$delimit     = NULL
inopts$check.names = F

rnames = inopts$rnames
inopts$rnames = NULL

inopts = c(inopts, list(row.names = if(rnames) 1 else NULL))
data   = read.table.inopts(infile, inopts)

eval(parse(text = {{args.helper | repr}}))

plot.bar(data, outfile, x, y, stacked, params, ggs, devpars = devpars)
