
library(optCluster)

d = read.table ({{in.infile | quote}}, sep="\t", header={{args.cnames | R}}, row.names={{args.rnames | lambda x: 'NULL' if not x else int(x)}}, check.names=F)
if ({{args.transpose | R}}) {
	d = t(d)
}
names = rownames(d)
{% if args.methods | lambda x: x == 'all' %}
#methods = c("agnes", "clara", "diana", "hierarchical", "kmeans", "pam", "som", "em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", "sa.poisson")
methods = c("agnes", "clara", "diana", "hierarchical", "pam", "em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", "sa.poisson")
{% elif args.methods | lambda x: isinstance(x, list) %}
methods = {{args.methods | Rvec}}
{% else %}
methods = unlist(strsplit({{args.methods | .replace(' ', '') | quote}}, ',', fixed = T))
{% endif %}

if ({{args.isCount | R}}) {
	methods = intersect(c("em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", "sa.poisson"), methods)
} else {
	methods = intersect(c("agnes", "clara", "diana", "hierarchical", "kmeans", "model", "pam", "som", "sota", "fanny"), methods)
}
minc = max({{args.minc}}, 2)
maxc = min({{args.maxc}}, nrow(d) - 1)

ret = optCluster(d, minc:maxc, clMethods = methods, countData = {{args.isCount | R}}, clVerbose = T, rankVerbose = F)


ret = optAssign(ret)

outfile = "{{out.outdir}}/clusters.txt"
metfile = "{{out.outdir}}/method.txt"
outfig  = "{{out.outdir}}/clusters.png"
outdir  = "{{out.outdir}}"
clusters = matrix(ret$cluster, ncol=1)
rownames(clusters) = names
write.table (clusters, outfile, row.names=T, sep="\t", quote=F, col.names=F)
write.table (ret$method, metfile, row.names=F, sep="\t", quote=F, col.names=F)
if ({{args.plot | R}}) {
	library(factoextra)
	pdata = list()
	pdata$data = d
	pdata$cluster = ret$cluster
	do.call(png, c(list(file = outfig), {{args.devpars | Rlist}}))
	print (fviz_cluster(pdata, data=d, main = ""))
	dev.off()
}
