
{{rimport}}('__init__.r', 'plot.r', 'sampleinfo.r')
library(methods)
options(stringsAsFactors = FALSE)

infile  = {{i.infile | quote}}
gfile   = {{i.gfile | quote}}
prefix  = {{i.infile | fn2 | quote}}
outdir  = {{o.outdir | quote}}
inopts  = {{args.inopts | R}}
tsform  = {{args.tsform or 'NULL'}}
filter  = {{args.filter or 'NULL'}}
plots   = {{args.plot | R}}
ggs     = {{args.ggs | R}}
params  = {{args.params | R}}
devpars = {{args.devpars | R}}

expr = read.table.inopts(infile, inopts)
if (is.function(filter)) {
	expr = filter(expr)
	outfile = file.path(outdir, basename(infile))
	write.table(expr, outfile, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
}
if (is.function(tsform)) {
	expr = tsform(expr)
}
saminfo = NULL
if (gfile != "") {
	saminfo = SampleInfo2$new(gfile)
	groups  = saminfo$all.groups()
}

if (plots$boxplot) {
	bpfile = file.path(outdir, paste0(prefix, '.boxplot.png'))
	plot.boxplot(expr, bpfile, stacked = F, devpars = devpars, params = params$boxplot, ggs = ggs$boxplot)
	if (!is.null(saminfo)) {
		for (group in groups) {
			bpfile = file.path(outdir, paste0(prefix, '.', group, '.boxplot.png'))
			plot.boxplot(
				expr[, saminfo$get.samples('Group', group), drop = FALSE], 
				bpfile, stacked = F, devpars = devpars, params = params$boxplot, ggs = ggs$boxplot)
		}
	}
}

if (plots$histogram) {
	histfile = file.path(outdir, paste0(prefix, ".histo.png"))
	plot.histo(stack(as.data.frame(expr)), histfile, devpars = devpars, params = params$histogram, ggs = ggs$histogram)
	if (!is.null(saminfo)) {
		for (group in groups) {
			histfile = file.path(outdir, paste0(prefix, '.', group, '.histo.png'))
			plot.histo(
				stack(expr[, saminfo$get.samples('Group', group), drop = FALSE]), 
				histfile, devpars = devpars, params = params$histogram, ggs = ggs$histogram)
		}
	}
}

if (plots$qqplot) {
	qqfile = file.path(outdir, paste0(prefix, ".qq.png"))
	plot.qq(expr, qqfile, stacked = FALSE, devpars = devpars, params = params$qqplot, ggs = ggs$qqplot)
	if (!is.null(saminfo)) {
		for (group in groups) {
			qqfile = file.path(outdir, paste0(prefix, '.', group, '.qq.png'))
			plot.qq(
				stack(expr[, saminfo$get.samples('Group', group), drop = FALSE]), 
				qqfile, stacked = FALSE, devpars = devpars, params = params$qqplot, ggs = ggs$qqplot)
		}
	}
}