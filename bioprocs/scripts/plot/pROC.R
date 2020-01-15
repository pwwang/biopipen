
{{rimport}}('__init__.r', 'plot.r')

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
outdir  = {{o.outdir | R}}
inopts  = {{args.inopts | R}}
params  = {{args.params | R}}
devpars = {{args.devpars | R}}
ggs     = {{args.ggs | R}}

indata = read.table.inopts(infile, inopts)

plotfile = file.path(outdir, '{{i.infile | stem}}.roc.png')
aucfile = file.path(outdir, '{{i.infile | stem}}.roc.txt')
params$returnTable = TRUE
ret = plot.roc(indata, plotfile, stacked = FALSE, params = params, ggs = ggs, devpars = devpars)
write.table(
	ret$table,
	aucfile,
	sep       = "\t",
	quote     = FALSE,
	col.names = TRUE,
	row.names = FALSE)
