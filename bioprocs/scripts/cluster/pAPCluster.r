
library ('factoextra')
library ('apcluster')
data = read.table ({{i.infile | quote}}, sep="\t", header={{args.cnames | R}}, row.names={{args.rnames | lambda x: 'NULL' if not x else int(x)}}, check.names=F)
if ({{args.transpose | R}}) {
	data = t(data)
}
data = data[, apply(data, 2, sd)!=0, drop=F]
ap        = apcluster(negDistMat(r=2), data)
clusters  = ap@clusters
nclust    = length (clusters)
clusts    = matrix (NA, ncol=1, nrow=ap@l)
rownames (clusts) = rownames(data)
for (i in 1:nclust) {
	cids  = unlist (clusters[i])
	for (cid in cids) clusts[cid, 1] = i
}
fvizobj   = {}
fvizobj$data = data
fvizobj$cluster = clusts

clustfile = file.path ({{o.outdir | quote}}, "clusters.txt")
write.table (clusts, clustfile, quote=F, col.names = F, sep="\t")
clustfig  = file.path ({{o.outdir | quote}}, "clusters.png")

do.call(png, c(list(file = clustfig), {{args.devpars | Rlist}}))
labels    = colnames(data)
print (fviz_cluster (fvizobj, main = ""))
dev.off()
