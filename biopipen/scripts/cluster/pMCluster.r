
library ('factoextra')
library ('mclust')
data = read.table ({{i.infile | quote}}, sep="\t", header={{args.cnames | R}}, row.names={{args.rnames | lambda x: 'NULL' if not x else int(x)}}, check.names=F)
if ({{args.transpose | R}}) {
	data = t(data)
}

data = data[, apply(data, 2, sd)!=0, drop=F]
mc        = Mclust (data, G={{args.minc}}:{{args.maxc}})
clustfile = file.path ({{o.outdir | quote}}, "clusters.txt")
write.table (mc$classification, clustfile, quote=F, col.names = F, sep="\t")
clustfig  = file.path ({{o.outdir | quote}}, "clusters.png")
do.call(png, c(list(file = clustfig), {{args.devpars | Rlist}}))
labels    = colnames(data)
print (fviz_cluster (mc, main = ""))
dev.off()
