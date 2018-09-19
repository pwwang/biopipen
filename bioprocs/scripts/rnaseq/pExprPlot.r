
{{rimport}}('__init__.r', 'plot.r')

infile  = {{i.infile | quote}}
prefix  = {{i.infile | fn2}}
outdir  = {{o.outdir | quote}}
transfm = {{args.transfm | lambda x: x or 'NULL'}}
plots   = {{args.plot | R}}
ggs     = {{args.ggs | R}}
devpars = {{args.devpars | R}}

expr = read.table.nodup(infile, header = T, row.names = 1, sep = "\t", check.names = F)

if (plot$boxplot) {
	bpfile = file.path(outdir, paste0(prefix, '.boxplot.png'))
	plot.boxplot(expr, bpfile, stack = T, devpars = devpars, ggs = ggs$boxplot)
}

if (plot$histogram) {
	histfile = file.path(outdir, paste0(prefix, ".histo.png"))
	plot.histo(stack(as.data.frame(expr)), histfile, devpars = devpars, ggs = ggs$histogram)
}
