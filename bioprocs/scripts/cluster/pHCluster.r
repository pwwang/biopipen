
data = read.table ({{i.infile | quote}}, sep="\t", header={{args.cnames | R}}, row.names={{args.rnames | lambda x: 'NULL' if not x else int(x)}}, check.names=F)
if ({{args.transpose | R}}) {
	data = t(data)
}
data = data[, apply(data, 2, sd)!=0, drop=F]

dmat = dist(data)
if ({{args.fast | R}}) {
	library('fastcluster')
}
orderfile = file.path ({{o.outdir | quote}}, "order.txt")
mergefile = file.path ({{o.outdir | quote}}, "merge.txt")
clustfig  = file.path ({{o.outdir | quote}}, "hclust.png")
hobj = hclust (dmat, method={{args.method | quote}})

do.call(png, c(list(file = clustfig), {{args.devpars | Rlist}}))
library('ggplot2')
library('ggdendro')
ggdendrogram (hobj, rotate={{args.rotate | R}})
dev.off()
orderobj = data.frame (Order=hobj$order, Labels=hobj$labels)
mergeobj = data.frame (Merge1=hobj$merge[,1], Merge2=hobj$merge[,2], Height=hobj$height)
write.table (orderobj, orderfile, quote=F, row.names=F, col.names = T, sep="\t")
write.table (mergeobj, mergefile, quote=F, row.names=F, col.names = T, sep="\t")

