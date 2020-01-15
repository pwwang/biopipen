{{rimport}}('plot.r', '__init__.r', 'sampleinfo.r')

library(methods)

set.seed(8525)

exprfile  = {{i.efile | R}}
groupfile = {{i.gfile | R}}
annofile  = {{args.annofile | R}}
cutoffway = {{args.pval | lambda x: str(x).split(':')[0] == 'q' and 'q' or 'p' | quote}}
cutoffval = {{args.pval | lambda x: str(x).split(':')[-1] | R}}
filters   = {{args.filter | lambda x: x if isinstance(x, (list, tuple)) else [int(y.strip()) for y in x.split(',')] | R}}
outdir    = {{o.outdir | R}}
allfile   = file.path(outdir, '{{i.efile | fn2}}-{{i.gfile | fn2}}.all.txt')
outfile   = {{o.outfile | R}}
tool      = {{args.tool | R}}
plot      = {{args.plot | R}}
ggs       = {{args.ggs | R}}
hmrows    = {{args.hmrows | R}}
devpars   = {{args.devpars | R}}

# get the exp data
ematrix   = read.table.nodup (exprfile,  header=T, row.names = 1, check.names=F, sep="\t")

# get groups
sampleinfo = SampleInfo2$new(groupfile)

gfactor   = factor(sampleinfo$all.groups())
group1    = levels(gfactor)[1]
group2    = levels(gfactor)[2]
samples1  = sampleinfo$get.samples(by = 'Group', value = group1)
samples2  = sampleinfo$get.samples(by = 'Group', value = group2)
samples   = c(samples1, samples2)
ematrix   = ematrix[, samples, drop=F]

n1 = length(samples1)
n2 = length(samples2)

ematrix = ematrix[rowSums(ematrix > filters[1]) >= filters[2], , drop=F]

# detect DEGs using limma
if (tool == 'limma') {
	# get model
	design = sampleinfo$as.edger.design()
	#print(design)
	library(limma)
	fit     = lmFit(ematrix, design)
	fit     = eBayes(fit)

	allgene = topTable(fit, n=nrow(ematrix), adjust="BH", coef = 2)
	write.table (pretty.numbers(allgene, list(P.Value..adj.P.Val = '%.2E', . = '%.3f')), allfile, quote=F, sep="\t")

	if (cutoffway == 'p') {
		degene  = allgene[allgene$P.Value < cutoffval,,drop=F]
	} else {
		degene  = allgene[allgene$adj.P.Val < cutoffval,,drop=F]
	}

	alllogFC = allgene$logFC
	allPval  = allgene$P.Value
	allQval  = allgene$adj.P.Val

	annodegs = NULL
	# add annotation
	if (annofile != '') {
		annos = read.table(annofile, header=F, sep="\t", row.names = 1, check.names = F )
		genes = data.frame(Gene = annos[rownames(degene), ])
		annodegs = cbind(genes, degene)
		write.table (pretty.numbers(annodegs, list(P.Value..adj.P.Val = '%.2E', logFC..AveExpr..t..B = '%.3f')), outfile, quote=F, sep="\t")
	} else {
		write.table (pretty.numbers(degene, list(P.Value..adj.P.Val = '%.2E', . = '%.3f')), outfile, quote=F, sep="\t")
	}

} else {
	stop(paste('Unsupported tool:', tool))
}

# MDS plot
if (plot$mdsplot) {
	library(limma)
	mdsplot = file.path (outdir, "mdsplot.png")
	do.call(png, c(list(filename=mdsplot), devpars))
	plotMDS(ematrix, col=c(rep("red2", n1), rep("blue2", n2)))
	dev.off()
}

# Volcano
if (plot$volplot != F) {
	volplot = file.path (outdir, "volcano.png")
	log2fc  = alllogFC
	fdr     = ifelse(cutoffway == 'p', allPval, allQval)
	plot.volplot(
		data.frame(logFC = log2fc, FDR = fdr),
		plotfile = volplot,
		fccut    = plot$volplot$fccut,
		pcut     = plot$volplot$pcut,
		ggs      = ggs$volplot,
		devpars  = devpars
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
		mat = ematrix[rownames(degene)[1:ngene], ]
		if (!is.null(annodegs)) {
			rownames(mat) = make.unique(as.character(as.vector(degene[1:ngene, 1])))
		}
		plot.heatmap(mat, hmap, ggs = ggs$heatmap, devpars = devpars)
	}
}

# maplots
if (plot$maplot) {
	maplot = file.path(outdir, "maplot.png")
	A = rowSums(ematrix) / ncol(ematrix)
	M = alllogFC
	threshold = ifelse(cutoffway == 'p', allPval < cutoffval, allQval < cutoffval)
	plot.maplot(
		data.frame(A, M),
		maplot,
		threshold,
		ggs = ggs$maplot,
		devpars = devpars
	)
}
