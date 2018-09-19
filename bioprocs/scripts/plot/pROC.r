
{{rimport}}('plot.r')

cnames  = as.logical({{args.cnames | R}})
rnames  = as.logical({{args.rnames | R}})
params  = {{args.params | R}}
devpars = {{args.devpars | R}}
ggs     = {{args.ggs | R}}
infile  = {{i.infile | R}}
outdir  = {{o.outdir | R}}
data = read.table(infile, sep = "\t", check.names = F, header = cnames, row.names = if (rnames) 1 else NULL)

if (!cnames) {
	colnames(data) = c('D', paste0('M', 1:(ncol(data)-1)))
}

params$returnAUC = T
plotfile = file.path(outdir, '{{i.infile | fn}}.roc.png')
aucfile  = file.path(outdir, '{{i.infile | fn}}.auc.txt')
aucs = plot.roc(data, plotfile, params, ggs, devpars)

write.table(aucs, aucfile, sep = "\t", quote = F)
