
{{rimport}}('__init__.r', 'plot.r')

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
outdir  = {{o.outdir | R}}
inopts  = {{args.inopts | R}}
params  = {{args.params | R}}
devpars = {{args.devpars | R}}
ggs     = {{args.ggs | R}}

indata = read.table.inopts(infile, inopts)

params$returnAUC = T
plotfile = file.path(outdir, '{{i.infile | fn}}.roc.png')
aucs = plot.roc(indata, plotfile, stacked = F, params = params, ggs = ggs, devpars = devpars)
aucs = t(as.data.frame(aucs))
write.table(pretty.numbers2(aucs, . = '%.3f'), outfile, sep = "\t", quote = F, col.names = F, row.names = T)
