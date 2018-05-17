{{rimport}}('plot.r', '__init__.r', 'sampleinfo.r')

library(methods)

exprfile  = {{in.efile | R}}
groupfile = {{in.gfile | R}}
pval      = {{args.pval | R}}
filters   = {{args.filter | lambda x: x if isinstance(x, (list, tuple)) else [int(y.strip()) for y in x.split(',')] | R}}
outdir    = {{out.outdir | R}}
allfile   = file.path(outdir, '{{in.efile | fn2}}-{{in.gfile | fn2}}.all.txt')
outfile   = {{out.outfile | R}}
tool      = {{args.tool | R}}
plot      = {{args.plot | R}}
ggs       = {{args.ggs | R}}
hmrows    = {{args.hmrows | R}}
devpars   = {{args.devpars | R}}

# get the exp data
ematrix   = read.table.nodup (exprfile,  header=T, row.names = 1, check.names=F, sep="\t")

# get groups
sampleinfo = SampleInfo(groupfile)$dataframe()

gfactor   = factor(sampleinfo$Group)
gfactor   = relevel(gfactor, as.vector(sampleinfo$Group[1]))
group1    = levels(gfactor)[1]
group2    = levels(gfactor)[2]
samples1  = row.names(sampleinfo[which(sampleinfo$Group == group1), 1, drop=F])
samples2  = row.names(sampleinfo[which(sampleinfo$Group == group2), 1, drop=F])
samples   = c(as.vector(samples1), as.vector(samples2))
ematrix   = ematrix[, samples, drop=F]

n1      = length(samples1)
n2      = length(samples2)

ematrix = ematrix[rowSums(ematrix > filters[1]) >= filters[2], , drop=F]
sicols  = colnames(sampleinfo)
sirows  = rownames(sampleinfo)
if ("Patient" %in% sicols && n1 != n2) {
	stop(paste0("Paired samples indicated, but number of samples is different in two groups (", n1, n2, ")."))
}

# detect DEGs using limma
if (tool == 'limma') {
	# get model
	if ("Patient" %in% sicols) {
		pairs  = factor(sampleinfo[which(sirows %in% samples), "Patient"])
		group  = factor(sampleinfo[which(sirows %in% samples), "Group"])
		group  = relevel(group, group2)
		design = model.matrix(~pairs + group)
	} else {
		group  = factor(sampleinfo[which(sirows %in% samples), "Group"])
		group  = relevel(group, group2)
		design = model.matrix(~group)
	}
	#print(design)
	library(limma)
	fit     = lmFit(ematrix, design)
	fit     = eBayes(fit)

	allgene = topTable(fit, n=nrow(ematrix), adjust="BH", coef = paste("group", group1, sep=""))
	write.table (allgene, allfile, quote=F, sep="\t")

	degene  = allgene[allgene$P.Value < pval,,drop=F]
	write.table (degene, outfile, quote=F, sep="\t")

	normedCounts = log2(ematrix + 1)
	alllogFC     = allgene$logFC
	allPval      = allgene$P.Value

} # leave it for futher tools

# MDS plot
if (plot$mdsplot) {
	library(limma)
	mdsplot = file.path (outdir, "mdsplot.png")
	do.call(png, c(list(filename=mdsplot), devpars))
	plotMDS(normedCounts, col=c(rep("red", n1), rep("blue", n2)))
	dev.off()
}

# Volcano
if (plot$volplot) {
	volplot = file.path (outdir, "volcano.png")
	log2fc  = alllogFC
	fdr     = allPval
	fdrcut  = rep(pval, length(allPval))
	plot.volplot(
		data.frame(logFC = log2fc, FDR = fdr, FDRCut = fdrcut),
		plotfile = volplot,
		ggs = ggs$volplot,
		devpars = devpars
	)
}

# heatmap
if (plot$heatmap) {
	if (nrow(degene) < 2) {
		logger('Cannot generate heatmap as < 2 DEGs detected.')
	} else {
		hmap   = file.path (outdir, "heatmap.png")
		ngene  = nrow(degene)
		if (hmrows <= 1) hmrows = as.integer(hmrows * ngene)

		ngene = min(hmrows, nrow(degene))
		mat   = normedCounts[rownames(degene[1:ngene, , drop=F]), ]
		plot.heatmap(mat, hmap, ggs = ggs$heatmap, devpars = devpars)
	}
}

# maplots
if (plot$maplot) {
	maplot = file.path(outdir, "maplot.png")
	A = rowSums(normedCounts) / ncol(normedCounts)
	M = alllogFC
	threshold = allPval < pval
	plot.maplot(
		data.frame(A, M, threshold),
		maplot,
		ggs = ggs$maplot,
		devpars = devpars
	)
}
