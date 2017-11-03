data     = read.table ("{{in.expfile}}", sep="\t", header={{args.header | R}}, row.names = NULL, check.names=F)
rnames   = make.unique(as.vector(data[,1]))
data[,1] = NULL
rownames(data)  = rnames

getGenelen = function() {
	gdata  = read.table({{args.refgene | quote}}, sep="\t", header=F, row.names=NULL, check.names=F)
	gnames = strsplit(as.vector(gdata[, ncol(gdata)]), ";", fixed=T)
	gnames = matrix(unlist(gnames), nrow=length(gnames[[1]]))[1,]
	gnames = strsplit(gnames, " ", fixed=T)
	gnames = matrix(unlist(gnames), nrow=length(gnames[[1]]))
	gnames = gnames[nrow(gnames),]
	gnames = make.unique(gnames)
	glen   = matrix(gdata[,5] - gdata[,4], ncol = 1)/1000
	colnames(glen) = c("GLEN")
	rownames(glen) = gnames
	return (glen)
}

{% if args.unit | lambda x: x == 'fpkm' or x == 'rpkm' %}
# convert to total 50,000,000 reads per sample
glen = getGenelen()
rn   = intersect(rownames(glen), rnames)
ret  = data[rn,,drop=F] * glen[rn,,drop=F]
ssum = colSums(ret)
ret  = 50000000 * ret / ssum

{% elif args.unit | lambda x: x == 'tpm' %}
glen = getGenelen()
rn   = intersect(rownames(glen), rnames)
ret  = data[rn,,drop=F]
ssum = colSums(ret)
ret  = 10000000 * ret / ssum
ret  = ret * glen[rn,,drop=F]
{% endif %}

write.table (ret, "{{out.outfile}}", quote=F, row.names=T, col.names={{args.header | R}}, sep="\t")

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
