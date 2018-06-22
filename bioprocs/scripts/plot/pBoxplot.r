{{rimport}}('plot.r')

library(ggplot2) # make aes available

infile  = {{in.infile | R}}
outfile = {{out.outfile | R}}
inopts  = {{args.inopts | R}}
devpars = {{args.devpars | R}}
x       = {{args.x | repr}}
y       = {{args.y | repr}}
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
data   = do.call(read.table, inopts)

if (!is.na(as.numeric(x))) x = as.integer(x)
if (!is.na(as.numeric(y))) y = as.integer(y)

eval(parse(text = {{args.helper | repr}}))
params  = {{args.params | R}}
ggs     = {{args.ggs | R}}

plot.boxplot(data, outfile, x, y, stacked, params, ggs, devpars = devpars)
