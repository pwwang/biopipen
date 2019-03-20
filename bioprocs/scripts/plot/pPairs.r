{{rimport}}('__init__.r', 'plot.r')
options(stringsAsFactors = TRUE)

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
inopts  = {{args.inopts | R }}
devpars = {{args.devpars | R}}
params  = {{args.params | R}}
ggs     = {{args.ggs | R}}

indata = read.table.inopts(infile, inopts)
plot.pairs(indata, outfile, params = params, ggs = ggs, devpars = devpars)
