
{{plot.boxplot.r}}
cnames = {{args.cnames | R}}
rnames = {{args.rnames | R}}
data = read.table ("{{in.datafile}}", sep="\t", header=cnames, row.names=NULL, check.names=F)
if (rnames) {
	rns   = make.unique(as.vector(data[,1]))
	data[,1] = NULL
	rownames(data)  = rns
}
plotBoxplot(data, {{out.outpng | quote}}, {{args.ggs | Rlist}}, {{args.devpars | Rlist}})
