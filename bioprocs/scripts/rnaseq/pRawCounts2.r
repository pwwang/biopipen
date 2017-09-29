data  = read.table ("{{in.expfile}}", sep="\t", header={{args.header | R}}, row.names = {{args.rownames}}, check.names=F)

if ("{{args.unit}}" == 'cpm') {
	library('edgeR')
	ret = cpm (data, log = {{args.log2 | R}})
} else if ("{{args.unit}}" == 'rpkm') {
	library('edgeR')
	genelen = read.table ("{{args.glenfile}}", header=F, row.names = 1, check.names = F)
	ret = rpkm (data, log = {{args.log2 | R}}, gene.length = as.vector(genelen))
} else {
	library('coseq')
	ret = transform_RNAseq(data, norm="TMM")
	ret = ret$normCounts
	if ({{args.log2 | R}})
		ret = log2(ret)
}

rnames = {{args.rownames}} == 1
cnames = {{args.header | R}}

write.table (ret, "{{out.outfile}}", quote=F, row.names=rnames, col.names=cnames, sep="\t")

# boxplot
{% if args.boxplot %}
{{ plotBoxplot }}
bpfile = file.path("{{out.outdir}}", "{{in.expfile | fn | fn}}.boxplot.png")
plotBoxplot(ret, bpfile, devpars = {{args.devpars | Rlist}}, ggs = {{args.boxplotggs | Rlist}})
{% endif %}

# heatmap
{% if args.heatmap %}
{{ plotHeatmap }}
hmfile = file.path("{{out.outdir}}", "{{in.expfile | fn | fn}}.heatmap.png")
hmexp  = if (nrow(ret) > {{args.heatmapn}}) ret[sample(nrow(ret),size={{args.heatmapn}}),] else ret
plotHeatmap(hmexp, hmfile, devpars = {{args.devpars | Rlist}}, ggs = {{args.heatmapggs | Rlist}})
{% endif %}

# histgram
{% if args.histplot %}
{{ plotHist }}
histfile = file.path("{{out.outdir}}", "{{in.expfile | fn | fn}}.hist.png")
plotHist(ret, histfile, devpars = {{args.devpars | Rlist}}, ggs = {{args.histplotggs | Rlist}})
{% endif %}
