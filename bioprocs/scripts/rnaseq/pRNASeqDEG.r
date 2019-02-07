{{rimport}}('plot.r', '__init__.r', 'sampleinfo.r')

library(methods)
options(stringsAsFactors = FALSE)
set.seed(8525)

exprfile  = {{i.efile | R}}
groupfile = {{i.gfile | R}}
cutoff    = {{args.cutoff | :{"by": "p", "value": a, "sign": "<"} if not isinstance(a, dict) else a | R }}
outdir    = {{o.outdir | R}}
allfile   = {{o.outfile | prefix | prefix | @append: ".all.txt" | quote}}
upfile    = {{o.outfile | prefix | prefix | @append: ".up.txt" | quote}}
downfile  = {{o.outfile | prefix | prefix | @append: ".down.txt" | quote}}
outfile   = {{o.outfile | R}}
tool      = {{args.tool | R}}
plots     = {{args.plot | R}}
ggs       = {{args.ggs | R}}
params    = {{args.params | R}}
devpars   = {{args.devpars | R}}
mapfile   = {{args.mapfile | quote}}

# get the exp data
ematrix   = read.table.nodup (exprfile,  header=T, row.names = 1, check.names=F, sep="\t")

# get groups
saminfo = SampleInfo2$new(groupfile, checkPaired = TRUE)

groups = saminfo$all.groups()
if (length(groups) != 2) {
	stop("I don't know how to do comparison between non-paired groups (# group must equal 2)")
}
group1   = groups[1]
group2   = groups[2]
samples1 = saminfo$get.samples(by = 'Group', value = group1)
samples2 = saminfo$get.samples(by = 'Group', value = group2)
# the real matrix
samples = c(samples1, samples2)
saminfo$select(samples)
ematrix = ematrix[, samples, drop = FALSE]

maps = NULL
if (mapfile != "") {
	maps = read.table(mapfile, header = FALSE, row.names = 1, sep = "\t", check.names = FALSE)
}

if (tool == 'edger') {
	# detect DEGs using edgeR
	# get model
	design = saminfo$as.edger.design()

	library(edgeR)
	dge = DGEList(counts = ematrix, group = saminfo$mat$Group)
	dge$samples$lib.size = colSums(dge$counts)
	dge = calcNormFactors(dge)

	disp = estimateDisp (dge, design)
	#fit  = glmQLFit (disp, design)
	#fit  = glmQLFTest (fit, coef = 2)
	fit = glmFit(disp, design)
	fit = glmLRT(fit, coef = 2)

	allgene = topTags (fit, n=nrow(fit$table), p.value = 1)
	allgene = allgene$table
	if (!is.null(maps)) {
		probes = rownames(allgene)
		allgene = cbind(Gene = maps[probes, 1], allgene)
		rownames(allgene) = probes
		rm(probes)
		colnames(allgene)[1] = 'Gene'
	}

	write.table (
		pretty.numbers2(
			allgene, PValue..FDR='%.2E', .='%.3f'
		), allfile, quote=FALSE, sep="\t")

	allgene = allgene[!is.na(allgene$PValue) & !is.na(allgene$FDR),, drop=F]
	if (cutoff$by == 'q') {
		degene = allgene[as.numeric(allgene$FDR) < cutoff$value,, drop = FALSE]
	} else {
		degene = allgene[as.numeric(allgene$PValue) < cutoff$value,, drop=FALSE]
	}
	write.table (
		pretty.numbers2(
			degene, PValue..FDR='%.2E', .='%.3f'
		), outfile, quote=FALSE, sep="\t")
	write.table (
		pretty.numbers2(
			degene[degene$logFC > 0,], PValue..FDR='%.2E', .='%.3f'
		), upfile, quote=FALSE, sep="\t")
	write.table (
		pretty.numbers2(
			degene[degene$logFC < 0,], PValue..FDR='%.2E', .='%.3f'
		), downfile, quote=FALSE, sep="\t")

	# for plots
	normedCounts = log2(dge$counts[rownames(allgene),, drop=FALSE] + 1)
	alllogFC     = allgene$logFC
	allPval      = ifelse(cutoff$by == 'q', allgene$FDR, allgene$PValue)

} else if (tool == 'deseq2') {
	# detect DEGs using DESeq2
	library(DESeq2)
	# model
	design = saminfo$as.deseq2.design()
	dge = DESeqDataSetFromMatrix(round(ematrix), design$coldata, design = design$design)

	dge     = DESeq(dge)
	allgene = results(dge)
	allgene = allgene[order(allgene$pvalue),,drop=F]
	if (!is.null(maps)) {
		probes = rownames(allgene)
		allgene = cbind(Gene = maps[probes, 1], allgene)
		rownames(allgene) = probes
		rm(probes)
		colnames(allgene)[1] = 'Gene'
	}

	write.table (
		pretty.numbers2(
			allgene, pvalue..padj='%.2E', .='%.3f'
		), allfile, quote=FALSE, sep="\t")

	allgene  = allgene[!is.na(allgene$pvalue) & !is.na(allgene$padj),,drop=F]
	if (cutoff$by == 'q') {
		degene = allgene[as.numeric(allgene$padj) < cutoff$value,, drop = FALSE]
	} else {
		degene = allgene[as.numeric(allgene$pvalue) < cutoff$value,, drop = FALSE]
	}
	write.table (
		pretty.numbers2(
			degene, pvalue..padj='%.2E', .='%.3f'
		), outfile, quote=FALSE, sep="\t")
	write.table (
		pretty.numbers2(
			degene[degene$log2FoldChange > 0, ], pvalue..padj='%.2E', .='%.3f'
		), upfile, quote=FALSE, sep="\t")
	write.table (
		pretty.numbers2(
			degene[degene$log2FoldChange < 0, ], pvalue..padj='%.2E', .='%.3f'
		), downfile, quote=FALSE, sep="\t")

	degene$logFC  = degene$log2FoldChange
	degene$PValue = degene$pvalue

	# for plots
	normedCounts  = log2(assay(dge)[rownames(allgene),, drop=FALSE] + 1)
	alllogFC      = allgene$log2FoldChange
	allPval       = ifelse(cutoff$by == 'q', allgene$padj, allgene$pvalue)
} else {
	stop(paste('Unsupported tool:', tool))
}

# MDS plot
if (plots$mdsplot) {
	library(limma)
	mdsplot = file.path (outdir, "mdsplot.png")
	do.call(png, c(list(filename=mdsplot), devpars))
	plotMDS(normedCounts, col=c(rep("red3", length(samples1)), rep("blue3", length(samples2))))
	dev.off()
}

# Volcano
if (plots$volplot) {
	volplot = file.path (outdir, "volcano.png")
	log2fc  = alllogFC
	fdr     = allPval
	plot.volplot(
		data.frame(logFC = log2fc, FDR = fdr),
		fccut    = params$volplot$fccut,
		pcut     = cutoff$value,
		plotfile = volplot,
		ggs      = ggs$volplot,
		devpars  = devpars
	)
}

# heatmap
if (plots$heatmap) {
	if (nrow(degene) < 2) {
		logger('Cannot generate heatmap as < 2 DEGs detected.')
	} else {
		hmap   = file.path (outdir, "heatmap.png")

		mat   = normedCounts[rownames(degene), ]
		plot.heatmap(mat, hmap, ggs = ggs$heatmap, devpars = devpars)
	}
}

# maplots
if (plots$maplot) {
	maplot = file.path(outdir, "maplot.png")
	A = rowSums(normedCounts) / ncol(normedCounts)
	M = alllogFC
	threshold = allPval < cutoff$value
	plot.maplot(
		data.frame(A, M),
		maplot,
		threshold,
		ggs = ggs$maplot,
		devpars = devpars
	)
}

# qqplot
if (plots$qqplot) {
	qqplot = file.path(outdir, 'qqplot.png')
	plot.qq(data.frame(`-log10(Pval)` = -log10(degene$PValue)), qqplot, stacked = TRUE, ggs = ggs$qqplot, devpars = devpars)
}
