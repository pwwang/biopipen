
{{rimport}}('plot.r', '__init__.r', 'sampleinfo.r')

exprfile  = {{i.expr      | R}}
batchfile = {{i.batch     | R}}
tool      = {{args.tool    | R}}
outfile   = {{o.outfile  | R}}
outdir    = {{o.outdir   | R}}
plot      = {{args.plot    | R}}
ggs       = {{args.ggs     | R}}
devpars   = {{args.devpars | R}}
hmrows    = {{args.hmrows  | R}}

expr    = read.table.nodup (exprfile, sep = "\t", header = TRUE, row.names = 1, check.names = F)

samples    = colnames(expr)
sampleinfo = SampleInfo(batchfile)
batch      = sampleinfo$dataframe(sampleinfo$select(sample = samples, get = c('Sample', 'Batch')))

batch      = factor(batch[samples, 'Batch'])

if (tool == 'combat') {
	library(sva)
	newexpr   = ComBat(dat=expr, batch=batch, par.prior = TRUE, mod = NULL)
	write.table (round(newexpr, 3), outfile, col.names=T, row.names=T, sep="\t", quote=F)
} else { # leave it for further tools
	stop('Unsupported tool: {{args.tool}}.')
}

# boxplot
if (plot$boxplot) {
	bpfile = file.path(outdir, "{{i.expr | fn2}}.boxplot.png")
	plot.boxplot(newexpr, bpfile, stack = T, devpars = devpars, ggs = ggs$boxplot)
}

# heatmap
if (plot$heatmap) {
	hmfile = file.path(outdir, "{{i.expr | fn2}}.heatmap.png")
	hmexp  = if (nrow(newexpr) > hmrows) newexpr[sample(nrow(newexpr),size=hmrows),] else newexpr
	plot.heatmap(hmexp, hmfile, devpars = devpars, ggs = ggs$heatmap)
}

# histgram
if (plot$histogram) {
	histfile = file.path(outdir, "{{i.expr | fn2}}.histo.png")
	plot.histo(as.data.frame(newexpr), histfile, devpars = devpars, ggs = ggs$histogram)
}
