{{rimport}}('__init__.r', 'plot.r')


infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
inopts  = {{args.inopts | R }}
devpars = {{args.devpars | R}}
x       = {{args.x | repr}}

data = read.table.inopts(infile, inopts)
if (!is.na(as.numeric(x))) x = as.integer(x)

eval(parse(text = {{args.helper | repr}}))
params  = {{args.params | R}}
ggs     = {{args.ggs | R}}

plot.density(data, outfile, x, params, ggs, devpars = devpars)
