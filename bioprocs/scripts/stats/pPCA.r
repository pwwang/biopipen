{{rimport}}('__init__.r', 'plot.r')

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
prefix   = {{o.outfile | prefix2 | quote}}
outdir   = {{o.outdir | quote}}
inopts   = {{args.inopts | R}}
annofile = {{args.anfile | quote}}
anopts   = {{args.anopts | R}}
devpars  = {{args.devpars | R}}
na       = {{args.na | R}}
seed     = {{args.seed | R}}
select   = {{args.select | R}}
set.seed(seed)

indata = read.table.inopts(infile, inopts)
if (!is.logical(na)) {
	indata[is.na(indata)] = na
}

indata = indata[complete.cases(indata), colSums(indata)!=0, drop = FALSE]

anno = NULL
if (annofile != "") {
	anno = read.table.inopts(annofile, anopts)[rownames(indata),,drop=FALSE]
}
plots = {{args.plots | R}}
pca   = prcomp(indata, scale = TRUE)

allpcfile    = paste0(prefix, '.allpcs.txt')
sdevfile     = paste0(prefix, '.sdev.txt')
centerfile   = paste0(prefix, '.center.txt')
rotationfile = paste0(prefix, '.rotation.txt')
write.table(round(pca$x, 3), allpcfile, sep = "\t", quote = FALSE)
sdev = data.frame(Sdev = pca$sdev, Percent = 100.0*pca$sdev/sum(pca$sdev))
sdev$CumPercent = cumsum(sdev$Percent)
rownames(sdev) = paste0('PC', 1:length(pca$sdev))
write.table(sdev, sdevfile, sep = "\t", quote = FALSE)
center = as.data.frame(pca$center)
colnames(center) = 'Center'
write.table(center, centerfile, sep = "\t", quote = FALSE)
rm(center)
write.table(pca$rotation, rotationfile, sep = "\t", quote = FALSE)

if (is.null(select)) {
	pcs = pca$x
} else if (select > 1) {
	pcs = pca$x[, 1:select, drop = FALSE]
} else {
	npcs = which(sdev$Percent >= select * 100)
	npcs = 1:max(c(2, npcs))
	pcs  = pca$x[, npcs, drop = FALSE]
}
write.table(round(pcs, 3), outfile, sep = "\t", quote = FALSE)

if (is.true(plots$scree)) {
	library(factoextra)
	screefile = paste0(prefix, '.screeplot.png')
	if (!is.list(plots$scree)) {
		plots$scree = list(ncp = 20)
	}
	do.call(png, c(list(screefile), devpars))
	print(do.call(fviz_screeplot, c(list(pca), plots$scree)))
	dev.off()
}

if (is.true(plots$var)) {
	library(factoextra)
	varfile = paste0(prefix, '.varplot.png')
	if (!is.list(plots$var)) {
		plots$var = list(repel = TRUE)
	}
	do.call(png, c(list(varfile), devpars))
	print(do.call(fviz_pca_var, c(list(pca), plots$var)))
	dev.off()
}

if (is.true(plots$bi)) {
	library(factoextra)
	bifile = paste0(prefix, '.biplot.png')
	if (!is.list(plots$bi)) {
		plots$bi = list(repel = TRUE)
	}
	do.call(png, c(list(bifile), devpars))
	print(do.call(fviz_pca_biplot, c(list(pca), plots$bi)))
	dev.off()
}

if (is.true(plots$clplot, 'any')) {
	library(factoextra)
	library(cluster)
	clusterfile = paste0(prefix, '.clusterplot.png')

	clmethod = plots$cluster$method
	plots$cluster$method = NULL
	plots$cluster$x = scale(pcs)
	cldata          = do.call(clmethod, plots$cluster)
	
	if (!is.list(plots$clplot)) {
		plots$clplot = list(repel = TRUE)
	}
	plots$clplot$object = cldata
	plots$clplot$data   = pcs
	plots$clplot$ggs    = NULL

	do.call(png, c(list(clusterfile), devpars))
	p = do.call(fviz_cluster, plots$clplot)
	if (!is.null(anno)) anno = cbind(p$data, anno)
	# reload anno
	ggs = {{args.plots.get('clplot', {}).get('ggs', {}) | R}}
	print(apply.ggs(p, ggs))
	dev.off()
}

if (is.true(plots$pairs)) {
	pcs = pca$x[, 1:plots$pairs$pcs, drop = FALSE]
	cname = NULL
	if (plots$pairs$anno == 'kmeans') {
		anno = as.data.frame(kmeans(pcs, plots$pairs$k, nstart = ifelse(is.null(plots$pairs$nstart), 25, plots$pairs$nstart), iter.max = ifelse(is.null(plots$pairs$niter), 1000, plots$pairs$niter))$clust)
		colnames(anno) = 'Cluster'
		pcs = cbind(pcs, anno)
		pcs$Cluster = as.factor(pcs$Cluster)
		cname = 'Cluster'
	} else if (file.exists(plots$pairs$anno)) {
		anno = read.table(plots$pairs$anno, row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)
		anno = anno[rownames(pcs),,drop = FALSE]
		anno[,1] = as.factor(anno[,1])
		cname = colnames(anno)[1]
		pcs = cbind(pcs, anno)
	}
	pairsfile = paste0(prefix, '.pairs.png')
	if (!is.null(cname)) {
		plot.pairs(pcs, pairsfile, params = list(mapping = ggplot2::aes_string(color = cname), upper = list(continuous = "density")), ggs = list(theme = list(axis.text.x = element_text(angle = 90, hjust = 1))))
	} else {
		plot.pairs(pcs, pairsfile, params = list(upper = list(continuous = "density")), ggs = list(theme = list(axis.text.x = element_text(angle = 90, hjust = 1))))
	}
}
