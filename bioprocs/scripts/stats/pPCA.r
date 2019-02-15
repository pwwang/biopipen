{{rimport}}('__init__.r', 'plot.r')

infile   = {{i.infile | quote}}
annofile = {{i.annofile | quote}}
outfile  = {{o.outfile | quote}}
prefix   = {{o.outfile | prefix2 | quote}}
outdir   = {{o.outdir | quote}}
inopts   = {{args.inopts | R}}
anopts   = {{args.anopts | R}}
devpars  = {{args.devpars | R}}
na       = {{args.na | R}}
seed     = {{args.seed | R}}
set.seed(seed)

indata = read.table.inopts(infile, inopts)
if (!is.logical(na)) {
	indata[is.na(indata)] = na
} else {
	indata = indata[complete.cases(indata),,drop = FALSE]
}
anno = NULL
if (annofile != "") {
	anno = read.table.inopts(annofile, anopts)[rownames(indata),,drop=FALSE]
}
plots = {{args.plots | R}}
pca   = prcomp(indata, scale = TRUE)

xfile        = outfile
sdevfile     = paste0(prefix, '.sdev.txt')
centerfile   = paste0(prefix, '.center.txt')
rotationfile = paste0(prefix, '.rotation.txt')
write.table(pca$x, xfile, sep = "\t", quote = FALSE)
sdev = data.frame(Sdev = pca$sdev, Percent = 100.0*pca$sdev/sum(pca$sdev))
sdev$CumPercent = cumsum(sdev$Percent)
rownames(sdev) = paste0('PC', 1:length(pca$sdev))
write.table(sdev, sdevfile, sep = "\t", quote = FALSE)
center = as.data.frame(pca$center)
colnames(center) = 'Center'
write.table(center, centerfile, sep = "\t", quote = FALSE)
rm(center)
write.table(pca$rotation, rotationfile, sep = "\t", quote = FALSE)

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

	clmethod             = plots$cluster$method
	npcs                 = plots$cluster$npcs
	plots$cluster$method = NULL
	plots$cluster$npcs   = NULL
	if (npcs > 1) {
		pcadata = pca$x[, 1:npcs, drop = FALSE]
	} else {
		pcs = which(sdev$CumPercent <= npcs * 100)
		if (length(pcs) < 2) pcs = 1:2
		pcadata = pca$x[, pcs, drop = FALSE]
	}
	plots$cluster$x = scale(pcadata)
	cldata          = do.call(clmethod, plots$cluster)
	
	if (!is.list(plots$clplot)) {
		plots$clplot = list(repel = TRUE)
	}
	plots$clplot$object = cldata
	plots$clplot$data   = pcadata
	plots$clplot$ggs    = NULL

	do.call(png, c(list(clusterfile), devpars))
	p = do.call(fviz_cluster, plots$clplot)
	if (!is.null(anno)) anno = cbind(p$data, anno)
	# reload anno
	ggs = {{args.plots.clplot.ggs | R}}
	print(apply.ggs(p, ggs))
	dev.off()
}
