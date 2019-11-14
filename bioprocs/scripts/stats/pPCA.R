{{rimport}}('__init__.r', 'plot.r')

set.seed(8525)
library(PCAtools)

infile   = {{i.infile | quote}}
metafile = {{i.metafile | quote}}
outfile  = {{o.outfile | quote}}
prefix   = {{o.outfile | prefix2 | quote}}
outdir   = {{o.outdir | quote}}
inopts   = {{args.inopts | R}}
devpars  = {{args.devpars | R}}
params   = {{args.params | R}}
na       = {{args.na | R}}
plots    = {{args.plots | R}}
npc      = {{args.npc | R}}

indata = read.table.inopts(infile, inopts)
if (!is.logical(na)) {
	indata[is.na(indata)] = na
}

indata = indata[complete.cases(indata),, drop = FALSE][rowSums(indata)!=0,,drop = FALSE]

metadata = NULL
if (metafile != "") {
	metadata = read.table.inopts(metafile, list(cnames = TRUE, rnames = TRUE))
}

metadata = metadata[colnames(indata),,drop = FALSE]

params$mat = indata
params$metadata = metadata
p = do.call(pca, params)

allpcfile   = paste0(prefix, '.allpcs.txt')
sdevfile    = paste0(prefix, '.sdev.txt')
rotatedfile = paste0(prefix, '.rotated.txt')

write.table(p$rotated, allpcfile, quote = FALSE, sep = "\t")
write.table(p$loadings, rotatedfile, quote = FALSE, sep = "\t")
write.table(p$sdev, sdevfile, quote = FALSE, sep = "\t")
sdev = data.frame(sdev = p$sdev, variance = p$variance)
rownames(sdev) = p$components
write.table(sdev, sdevfile, quote = FALSE, sep = "\t")

# try select optimal number of PCs and save them into output file
npcmethod = NULL
if (npc == 'horn') {
	npcmethod = npc
	horn <- parallelPCA(indata)
	npc = horn$n
} else if (npc == 'elbow') {
	npcmethod = npc
	npc = findElbowPoint(p$variance)
} else if (!is.numeric(npc)) {
	stop('Expected a number, "horn" or "elbow" for optimal number of PCs.')
}
write.table(p$rotated[, 1:npc], outfile, quote = FALSE, sep = "\t")

# start plotting
if (is.list(plots$scree) || plots$scree != FALSE) {
	screepng = paste0(prefix, '.scree.png')
	if (!is.list(plots$scree)) {
		plots$scree = list()
	}
	do.call(png, c(list(filename = screepng), devpars))
	screeggs = plots$scree$ggs
	plots$scree$ggs = NULL
	plots$scree$pcaobj = p
	plots$scree$components = getComponents(p, 1:20)
	plots$scree$vline = npc
	if (!is.null(npcmethod)) {
		screeggs = c(screeggs, list(
			geom_text = list(aes(npc+1, 50, label = npcmethod, vjust = -1))
		))
	}
	print(apply.ggs(do.call(screeplot, plots$scree), screeggs))
	dev.off()
}

if (is.list(plots$bi) || plots$bi != FALSE) {
	bipng = paste0(prefix, '.bi.png')
	if (!is.list(plots$bi)) {
		plots$bi = list()
	}
	do.call(png, c(list(filename = bipng), devpars))
	plots$bi$pcaobj = p
	print(do.call(biplot, plots$bi))
	dev.off()
}

if (is.list(plots$pairs) || plots$pairs != FALSE) {
	pairspng = paste0(prefix, '.pairs.png')
	if (!is.list(plots$pairs)) {
		plots$pairs = list()
	}
	do.call(png, c(list(filename = pairspng), devpars))
	plots$pairs$pcaobj = p
	plots$pairs$components = getComponents(p, 1:5)
	print(do.call(pairsplot, plots$pairs))
	dev.off()
}

if (is.list(plots$loadings) || plots$loadings != FALSE) {
	loadingspng = paste0(prefix, '.loadings.png')
	if (!is.list(plots$loadings)) {
		plots$loadings = list()
	}
	do.call(png, c(list(filename = loadingspng), devpars))
	plots$loadings$pcaobj = p
	print(do.call(plotloadings, plots$loadings))
	dev.off()
}

if (is.list(plots$eigencor) || plots$eigencor != FALSE) {
	eigencorpng = paste0(prefix, '.eigencor.png')
	if (!is.list(plots$eigencor)) {
		plots$eigencor = list()
	}
	do.call(png, c(list(filename = eigencor), devpars))
	plots$eigencor$pcaobj = p
	print(do.call(eigencorplot, plots$eigencor))
	dev.off()
}
