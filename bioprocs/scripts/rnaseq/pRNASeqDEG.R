{{rimport}}('plot.r', '__init__.r', 'sampleinfo.r')

library(methods)
options(stringsAsFactors = FALSE)
set.seed(8525)

exprfile  = {{i.efile | R}}
groupfile = {{i.gfile | R}}
cutoff    = {{args.cutoff | ?isinstance: dict
						  | :{"by": _.get("by", "p"), "value": _.get("value", .05)}
						  | :{"by": "p", "value": _, "sign": "<"} | R }}
outdir    = {{o.outdir | R}}
allfile   = {{o.outfile | prefix | prefix | @append: ".all.xls" | quote}}
upfile    = {{o.outfile | prefix | prefix | @append: ".up.xls" | quote}}
downfile  = {{o.outfile | prefix | prefix | @append: ".down.xls" | quote}}
outfile   = {{o.outfile | R}}
tool      = {{args.tool | R}}
inopts    = {{args.inopts | R}}
ggs       = {{args.ggs | R}}
params    = {{args.params | R}}
devpars   = {{args.devpars | R}}
mapping   = {{args.mapping | quote}}

# get the exp data
indata = read.table.inopts(exprfile, inopts)

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
ematrix = indata[, samples, drop = FALSE]

maps = NULL
replacecol = NULL
if (is.true(mapping)) {
	splits = unlist(strsplit(mapping, ':'))
	if (length(splits) == 1) {
		if (file.exists(splits)) {
			maps = read.table.inopts(mapping, list(rnames = TRUE, colnames = TRUE))
		} else {
			replacecol = ifelse(is.na(as.numeric(splits)), splits, as.numeric(splits))
		}
	} else {
		replacecol = ifelse(is.na(as.numeric(splits[2])), splits[2], as.numeric(splits[2]))
		maps = read.table.inopts(mapping, list(rnames = TRUE, colnames = is.na(as.numeric(replacecol))))
	}
} else if (is.null(maps) && ncol(indata) > length(samples)) {
	maps = indata[, setdiff(colnames(indata), samples), drop = FALSE]
}

do_edger = function() {
	require('edgeR')
	# detect DEGs using edgeR
	# get model
	design = saminfo$as.edger.design()
	dge = DGEList(counts = ematrix, group = saminfo$mat$Group)
	dge$samples$lib.size = colSums(dge$counts)
	dge = calcNormFactors(dge)
	disp = estimateDisp (dge, design)
	#fit  = glmQLFit (disp, design)
	#fit  = glmQLFTest (fit, coef = 2)
	fit = glmFit(disp, design)
	fit = glmLRT(fit, coef = 2)

	allgene  = topTags(fit, n=nrow(fit$table), p.value = 1)
	allgenes = rownames(allgene)
	allgene  = data.frame(
		logFC  = allgene$table$logFC,
		PValue = allgene$table$PValue,
		FDR    = allgene$table$FDR,
		logCPM = allgene$table$logCPM,
		LR = allgene$table$LR
	)
	rownames(allgene) = allgenes
	allgene = allgene[!is.na(allgene$PValue) & !is.na(allgene$FDR),,drop=FALSE]
	list(allgene = allgene, normcounts = log2(dge$counts[rownames(allgene),, drop=FALSE] + 1))
}

do_deseq2 = function() {
	# detect DEGs using DESeq2
	library(DESeq2)
	# model
	design = saminfo$as.deseq2.design()
	dge = DESeqDataSetFromMatrix(round(ematrix), design$coldata, design = design$design)

	dge      = DESeq(dge)
	allgene  = results(dge)
	allgene  = as.data.frame(allgene[order(allgene$pvalue),,drop=FALSE])
	allgenes = rownames(allgene)
	allgene  = data.frame(
		logFC    = allgene$log2FoldChange,
		PValue   = allgene$pvalue,
		FDR      = allgene$padj,
		baseMean = allgene$baseMean,
		lfcSE    = allgene$lfcSE,
		stat     = allgene$stat
	)
	rownames(allgene) = allgenes
	allgene = allgene[!is.na(allgene$PValue) & !is.na(allgene$FDR),,drop=FALSE]
	list(allgene = allgene, normcounts = log2(assay(dge)[rownames(allgene),, drop=FALSE] + 1))
}

if (tool == 'edger') {
	ret = do_edger()
} else if (tool == 'deseq2') {
	ret = do_deseq2()
} else {
	stop(paste('Unsupported tool:', tool))
}

# save allgenes
# add maps to allgenes
if (!is.null(maps)) {
	probes = rownames(ret$allgene)
	if (is.true(replacecol)) {
		rownames(ret$allgene) = maps[, replacecol]
		restmaps = ifelse(is.numeric(replacecol),
			colnames(maps)[-replacecol],
			setdiff(colnames(maps), replacecol))
		allgene = cbind(ret$allgene, probes, maps[, restmaps, drop = FALSE])
		colnames(allgene)[ncol(allgene) + 1] = 'Original_Probe'
	} else {
		allgene = cbind(ret$allgene, maps[probes, , drop=FALSE])
	}
}
write.xls(allgene, allfile)

# get and save degs
if (cutoff$by == 'q') {
	degene = allgene[allgene$FDR < cutoff$value,,drop = FALSE]
} else {
	degene = allgene[allgene$PValue < cutoff$value,,drop = FALSE]
}
write.xls(degene, outfile)

# save up genes
write.xls(degene[degene$logFC>0,,drop=FALSE], upfile)
# save down genes
write.xls(degene[degene$logFC<0,,drop=FALSE], downfile)

# MDS plot
if (!is.logical(ggs$mdsplot) || ggs$mdsplot!=FALSE) {
	mdsplot = file.path (outdir, "mdsplot.png")
	mdsdata = data.frame(Group = c(rep(group1, length(samples1)), rep(group2, length(samples2))))
	rownames(mdsdata) = c(samples1, samples2)
	mdsdata = cbind(mdsdata, t(ret$normcounts[rownames(degene), c(samples1, samples2),drop=FALSE]))
	plot.mds(mdsdata, mdsplot, ggs = ggs$mdsplot, devpars = devpars)
}

# Volcano
if (!is.logical(ggs$volplot) || ggs$volplot!=FALSE) {
	volplot = file.path (outdir, "volcano.png")
	plot.volplot(
		allgene[, c('logFC', 'FDR')],
		volplot,
		params = list(logfccut = list.get(params$volplot, 'logfccut', 2),
					  pcut = list.get(params$volplot, 'pcut', 0.05),
					  hilight = list.get(params$volplot, 'hilight', NULL)),
		ggs      = ggs$volplot,
		devpars  = devpars
	)
}

# heatmap
if (!is.logical(ggs$heatmap) || ggs$heatmap!=FALSE) {
	if (nrow(degene) < 2) {
		logger('Cannot generate heatmap as < 2 DEGs detected.')
	} else {
		hmap   = file.path (outdir, "heatmap.png")
		ngenes = list.get(params$heatmap, 'ngenes', 500, check.names = TRUE)
		hmgenes = rownames(degene)
		if (!is.null(ngenes) && nrow(degene) > ngenes*2) {
			logger('Too many degs for heatmap, will cut it to', ngenes * 2)
			hmgenes = c(rownames(degene[degene$logFC > 0,,drop=FALSE])[1:ngenes],
						rownames(degene[degene$logFC < 0,,drop=FALSE])[1:ngenes])
		}
		mat   = ret$normcounts[hmgenes, ]
		params$heatmap$ngenes = NULL
		plot.heatmap2(mat, hmap, params = params$heatmap, draw = ggs$heatmap, devpars = devpars)
	}
}

# maplots
if (!is.logical(ggs$maplot) || ggs$maplot!=FALSE) {
	maplot = file.path(outdir, "maplot.png")
	A = rowSums(ret$normcounts) / ncol(ret$normcounts)
	M = allgene$logFC
	threshold = allgene$PValue < cutoff$value
	plot.maplot(
		data.frame(A, M),
		maplot,
		threshold,
		ggs = ggs$maplot,
		devpars = devpars
	)
}

# qqplot
if (!is.logical(ggs$qqplot) || ggs$qqplot!=FALSE) {
	qqplot = file.path(outdir, 'qqplot.png')
	plot.qq(data.frame(`-log10(Pval)` = -log10(degene$PValue)), qqplot, ggs = ggs$qqplot, devpars = devpars)
}
