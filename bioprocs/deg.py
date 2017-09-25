from pyppl import Proc
from .utils import cbindFill, plot

"""
@name:
	pExpdirMatrix
@description:
	Convert expression files to expression matrix
	File names will be used as sample names (colnames)
	Each gene and its expression per line.
@input:
	`expdir:file`:  the directory containing the expression files, could be gzipped
@output:
	`expfile:file`: the expression matrix
"""
pExpdir2Matrix = Proc()
pExpdir2Matrix.input           = "expdir:file"
pExpdir2Matrix.output          = "expfile:file:{{expdir | fn}}.matrix.txt"
pExpdir2Matrix.lang            = "Rscript"
pExpdir2Matrix.args.pattern    = '*'
pExpdir2Matrix.args.header     = False
pExpdir2Matrix.args.excl       = ["^Sample", "^Composite", "^__"]
pExpdir2Matrix.tplenvs.cbindFill = cbindFill.r
pExpdir2Matrix.script          = """
setwd("{{in.expdir}}")

{{cbindFill}}

isGoodRname = function(rname) {
	for (excl in {{args.excl | Rvec}}) {
		if (grepl(excl, rname)) {
			return (FALSE)
		}
	}
	return (TRUE)
}
isGoodRname = Vectorize(isGoodRname)

exp = c()
for (efile in Sys.glob({{args.pattern | quote}})) {
	write(paste("pyppl.log: Reading", efile, "..."), stderr())
	sample = tools::file_path_sans_ext(basename(efile))
	if (grepl ('.gz$', efile)) efile = gzfile (efile)
	tmp    = read.table (efile, sep="\\t", header=F, row.names = 1, check.names=F)
	rnames = rownames(tmp)
	rnames = rnames[isGoodRname(rnames)]
	tmp    = tmp[rnames,,drop=F]
	colnames (tmp) = c(sample)
	exp    = cbindFill (exp, tmp)
}

write.table (exp, "{{out.expfile}}", col.names=T, row.names=T, sep="\\t", quote=F)
"""

"""
	`bcvplot`: whether to plot biological coefficient of variation, default: True
	`displot`: whether to plot biological coefficient of variation, default: True
	`fcplot`:  whether to plot fold changes, default: True
"""
pRseqDEG              = Proc(desc = 'Detect DEGs by RNA-seq data.')
pRseqDEG.input        = "efile:file, group1, group2"
pRseqDEG.output       = [
	"outfile:file:{{in.group1 | .split(':', 1)[0]}}-{{in.group2 | .split(':', 1)[0]}}.degs/{{in.group1 | .split(':', 1)[0]}}-{{in.group2 | .split(':', 1)[0]}}.degs.txt", 
	"outdir:dir:{{in.group1 | .split(':', 1)[0]}}-{{in.group2 | .split(':', 1)[0]}}.degs"
]
pRseqDEG.args.tool    = 'edger' # limma, deseq2
pRseqDEG.args.filter  = '1,2'
pRseqDEG.args.pval    = 0.05
pRseqDEG.args.paired  = False
pRseqDEG.args.bcvplot = True
pRseqDEG.args.displot = True
pRseqDEG.args.fcplot  = True
pRseqDEG.args.heatmap = False
pRseqDEG.args.listall = False
pRseqDEG.args.batches = False # True to be the same as the groups, otherwise, list of batches: ["s1,s2,...", "s3,s4,..."]
pRseqDEG.tplenvs.plotHeatmap = plot.heatmap.r
pRseqDEG.lang         = "Rscript"
pRseqDEG.script       = """
# require sva to remove batch effect
require('sva')
trim <- function (x) gsub("^\\\\s+|\\\\s+$", "", x)
# get the exp data
ematrix    = read.table ("{{in.efile}}",  header=T, row.names = 1, check.names=F, sep="\\t")
cnames     = colnames(ematrix)

# get group names
groups1    = unlist(strsplit("{{in.group1}}", ":", fixed = TRUE))
groups2    = unlist(strsplit("{{in.group2}}", ":", fixed = TRUE))
group1name = groups1[1]
group2name = groups2[1]
group1     = unlist(strsplit(groups1[2], ",", fixed = TRUE))
group2     = unlist(strsplit(groups2[2], ",", fixed = TRUE))
group1idx  = as.integer (group1)
group2idx  = as.integer (group2)
if (NA %in% group1idx) {
	group1idx = match(group1, cnames)
} 
if (NA %in% group2idx) {
	group2idx = match(group2, cnames)
} 

ematrix    = ematrix[, c(group1idx, group2idx), drop=FALSE]
group1idx  = 1:length(group1idx)
group2idx  = (length(group1idx) + 1):ncol(ematrix)
cnames     = colnames(ematrix)
n1         = length(group1idx)
n2         = length(group2idx)
paired     = {{args.paired | R}}
pval       = {{args.pval | lambda x: float(x)}}
filters    = c({{args.filter}})
tool       = {{args.tool | quote}}

if (paired && n1 != n2) {
	stop(paste("Flag paired is TRUE, but number of samples is different in two groups (", n1, n2, ").", sep=""))
}

# get batches
batches    = NULL
argBatches = {{args.batches | lambda x: 'TRUE' if x is True else 'FALSE' if x is False else 'c('+ ','.join(['"'+ v +'"' for v in x]) +')'}}
print(paste('pyppl.log:', argBatches), stderr())
if (!is.logical(argBatches)) {
	batches = vector(length = length(cnames))
	for (i in 1:length(argBatches)) {
		b     = argBatches[i]
		bname = paste('B', i, sep='')
		ss    = sapply(unlist(strsplit(b, ",", fixed=T)), trimws)
		sidx  = as.integer(ss)
		if (NA %in% sidx) {
			sidx = match(ss, cnames)
		}
		batches[sidx] = bname
	}
} else if (argBatches) {
	batches = c(rep(group1name, n1), rep(group2name, n2))
}

pairs = vector(mode="numeric")
group = vector(mode="character")

if (tool == "edger") {
	library('edgeR')
	if (paired) {
		for (i in 1:n1) {
			pairs = c(pairs, i, i)
			group = c(group1name, group2name)
		}
		pairs  = factor(pairs)
		group  = factor(group)
		design = model.matrix(~pairs + group)
	} else {
		group  = c(rep(group1name, n1), rep(group2name, n2))
		group  = factor(group)
		design = model.matrix(~group)
	}
	
	dge    = DGEList(counts = ematrix, group = group)
	dge    = dge[rowSums(cpm(dge)>filters[1]) >= filters[2], ]
	dge$samples$lib.size = colSums(dge$counts)
	dge    = calcNormFactors(dge)
	
	if (!is.null(batches)) {
		mod1 = model.matrix(~batches)
		mod0 = mod1[,1]
		svas   = svaseq(cpm(dge), mod1, mod0)
		design = cbind(design, svas$sv)
	}
	
	disp   = estimateDisp (dge, design)
	fit    = glmFit (disp, design)
	fit    = glmLRT (fit)
	
	degs   = topTags (fit, n=nrow(fit$table), p.value = {{args.pval}})
	write.table (degs$table, "{{outfile}}", quote=F, sep="\\t")
	degs$names = rownames(degs$table)
	
	if ({{args.listall | R}}) {
		allgenes = topTags (fit, n=nrow(fit$table), p.value = 1)
		write.table (allgenes$table, "{{out.outdir}}/{{outfile | fn | fn}}.all.txt", quote=F, sep="\\t")
	}
	
	if ({{args.bcvplot | R}}) {
		bcvplot = file.path ("{{out.outdir}}", "bcvplot.png")
		png (file=bcvplot)
		plotMDS (dge, method="bcv", col=as.numeric(dge$samples$group))
		legend("bottomleft", as.character(unique(dge$samples$group)), col=2:1, pch=20)
		dev.off()
	}

	if ({{args.displot | R}}) {
		displot = file.path ("{{out.outdir}}", "displot.png")
		png (file=displot)
		plotBCV (disp)
		dev.off()
	}

	if ({{args.fcplot | R}}) {
		deg    = decideTestsDGE(fit, p.value = {{args.pval}})
		fcplot = file.path ("{{out.outdir}}", "fcplot.png")
		png (file=fcplot)
		tags = rownames(disp)[as.logical(deg)]
		plotSmear (fit, de.tags=tags)
		abline(h = c(-2, 2), col = "blue")
		dev.off()
	}
}

if ({{args.heatmap | R}}) {
	hmap = file.path ("{{out.outdir}}", "heatmap.png")
	{{plotHeatmap}}
	tmatrix    = ematrix[degs$names, group1idx]
	tmatrix    = tmatrix + 1
	nmatrix    = ematrix[degs$names, group2idx]
	nmatrix    = rowSums(nmatrix)/ncol(nmatrix)
	nmatrix    = nmatrix + 1
	log2fc     = log2(apply(tmatrix, 2, function(col) col/nmatrix))
	plotHeatmap(log2fc, list(filename=hmap, show_rownames=FALSE))
}


"""

"""
@name:
	pDEGByEdgeR
@description:
	Call DEG from expression matrix
@input:
	`expfile:file`: the expression matrix
	`group1`:       columns of group1 (separated by comma)
	`group2`:       columns of group2 (separated by comma)
	`group1name`:   the name of group1
	`group2name`:   the name of group2   
@output:
	`degdir:dir`:   the output directory containing DEGs and plots
@args:
	`filter`:  the pair (X,Y) on how to filter the data by cpm (`d <- d[rowSums(cpm(d)>X) >= Y,]`). Default: "1,2"
		- keep genes with at least X counts per million (CPM) in at least Y samples
	`pval`:    the cutoff of DEGs (default: .05)
	`paired`:  whether the samples are paired, default: False
	`bcvplot`: whether to plot biological coefficient of variation, default: True
	`displot`: whether to plot biological coefficient of variation, default: True
	`fcplot`:  whether to plot fold changes, default: True
@requires:
	[edgeR](https://bioconductor.org/packages/release/bioc/html/edger.html)
"""
pDEGByEdgeR = Proc()
pDEGByEdgeR.input     = "expfile:file, group1, group2, group1name, group2name"
pDEGByEdgeR.output    = "degdir:dir:{{expfile | fn}}.{{in.group1name}}-{{in.group2name}}.deg"
pDEGByEdgeR.args      = {'filter': "1,2", 'pval': 0.05, 'paired': False, 'bcvplot': True, 'displot': True, 'fcplot': True}
pDEGByEdgeR.defaultSh = "Rscript"
pDEGByEdgeR.script    = """
library('edgeR')
group1  = unlist(strsplit("{{in.group1}}", ","))
group2  = unlist(strsplit("{{in.group2}}", ","))

expmatrix = read.table ("{{expfile}}",  header=T, row.names = 1, check.names=F, sep="\\t")
g1names   = match (group1, colnames(expmatrix))
g2names   = match (group2, colnames(expmatrix))
expmatrix = expmatrix [, c(g1names, g2names)]
pairs = vector(mode="numeric")
group = vector(mode="character")
ng1   = length(g1names)
ng2   = length(g2names)
if ({{args.paired | R}}) {
	for (i in 1:ng1) {
		pairs = c(pairs, i, i)
		group = c(group, "{{in.group1name}}", "{{in.group2name}}")
	}
	pairs  = factor(pairs)
	group  = factor(group)
	design = model.matrix(~pairs+group)
} else {
	group  = c(rep("{{in.group1name}}", ng1), rep("{{in.group2name}}", ng2))
	group  = factor(group)
	design = model.matrix(~group)
}

# filter
dobj   = DGEList(counts=expmatrix, group=group)
filter = noquote(unlist(strsplit("{{args.filter}}", ",")))
fX     = as.numeric (filter[1])
fY     = as.numeric (filter[2])
dobj   = dobj[rowSums(cpm(dobj)>fX) >= fY, ]
dobj$samples$lib.size = colSums(dobj$counts)

# normalize
dobj = calcNormFactors(dobj, method="TMM")

if ({{args.bcvplot | R}}) {
	bcvplot = file.path ("{{degdir}}", "bcvplot.png")
	png (file=bcvplot)
	plotMDS (dobj, method="bcv", col=as.numeric(dobj$samples$group))
	legend("bottomleft", as.character(unique(dobj$samples$group)), col=2:1, pch=20)
	dev.off()
}

disp <- estimateDisp (dobj, design)
if ({{args.displot | R}}) {
	displot = file.path ("{{degdir}}", "displot.png")
	png (file=displot)
	plotBCV (disp)
	dev.off()
}

fit    = glmFit (disp, design)
fit    = glmLRT (fit)

if ({{args.fcplot | R}}) {
	deg    = decideTestsDGE(fit, p.value = {{args.pval}})
	fcplot = file.path ("{{degdir}}", "fcplot.png")
	png (file=fcplot)
	tags = rownames(disp)[as.logical(deg)]
	plotSmear (fit, de.tags=tags)
	abline(h = c(-2, 2), col = "blue")
	dev.off()
}

out    = topTags (fit, n=nrow(fit$table), p.value = {{args.pval}})
write.table (out$table, file.path("{{degdir}}", "degs.txt"), quote=F, sep="\\t")
"""

"""
@name:
	pMArrayLimma
@description:
	Call degs of microarray data by limma
@input:
	`expfile:file`: The expression matrix
	`group1`:      The 1st group
	`group2`:      The 2nd group
	`group1name`:  The name of 1st group
	`group2name`:  The name of 2nd group
@output:
	`degdir:dir`:   the output directory containing DEGs and plots
@args:
	`norm`:    the normalization methods, separated by comma. Support normalization methods: `quan` (quantile). Default: "quan"
	`boxplot`: draw boxplot? Default: True
	`paired`:  whether the samples are paired, default: False
	`filter`:  the pair (X,Y) on how to filter the data by expression (`d <- d[rowSums(d>X) >= Y,]`). Default: "1,2"
		- keep genes with at least X exp in at least Y samples
	`pval`:    the pvalue cutoff of DEGs (default: .05)
	`qval`:    the qvalue cutoff of DEGs (default: .05)
	`heatmap`: whether to plot heatmap, default: True
	`hmn`:     Number of gene used for heatmap, default: 50
	`hmmar`:   Margins for heatmap, default: "10,7"
	`volplot`: whether to plot the volcano plot, default: True
@requires:
	[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
	[ggplot2](https://bioconductor.org/packages/release/bioc/html/ggplot2.html)
"""
pMArrayLimma = Proc()
pMArrayLimma.input  = "expfile:file, group1, group2, group1name, group2name"
pMArrayLimma.output = "degdir:dir:{{expfile | fn}}.{{in.group1name}}-{{in.group2name}}.deg"
pMArrayLimma.args   = {'norm': "quan", 'boxplot': True, 'paired': False, "filter": "1,2", "qval": .05, "pval": .05, "heatmap": True, "hmn": 50, "volplot": True, "hmmar": "10,7"}
pMArrayLimma.lang   = "Rscript"
pMArrayLimma.script = """
library(limma)
library(ggplot2)
group1  = unlist(strsplit("{{in.group1}}", ","))
group2  = unlist(strsplit("{{in.group2}}", ","))

expmatrix = read.table ("{{expfile}}",  header=T, row.names = 1, check.names=F, sep="\\t")
g1names   = match (group1, colnames(expmatrix))
g2names   = match (group2, colnames(expmatrix))
expmatrix = expmatrix [, c(g1names, g2names)]

# filter
filter     = noquote(unlist(strsplit("{{args.filter}}", ",")))
fX         = as.numeric (filter[1])
fY         = as.numeric (filter[2])
expmatrix  = expmatrix[rowSums(expmatrix>fX) >= fY, ]

if ("quan" %in% unlist(strsplit("{{args.norm}}", ","))) {
	expmatrix = normalizeBetweenArrays (expmatrix)
}

if ({{args.boxplot | R}}) {
	png (file = file.path("{{degdir}}", "{{expfile | fn}}.{{in.group1name}}-{{in.group2name}}.boxplot.png"), res=300, width=2000, height=2000)
	boxplot (expmatrix, las=2)
	dev.off()	
}

pairs = vector(mode="numeric")
group = vector(mode="character")
ng1   = length(g1names)
ng2   = length(g2names)
if ({{args.paired | R}}) {
	for (i in 1:ng1) {
		pairs = c(pairs, i, i)
		group = c(group, "{{in.group1name}}", "{{in.group2name}}")
	}
	pairs  = factor(pairs)
	group  = factor(group)
	design = model.matrix(~pairs+group)
} else {
	group  = c(rep("{{in.group1name}}", ng1), rep("{{in.group2name}}", ng2))
	group  = factor(group)
	design = model.matrix(~group)
}

fit    = lmFit (expmatrix, design)
fit    = eBayes(fit)
out    = topTable (fit, number = nrow(expmatrix), sort.by="p")
out    = out[out$adj.P.Val<{{args.qval}} & out$P.Value<{{args.pval}},]

degfile = file.path("{{degdir}}", "{{expfile | fn}}.{{in.group1name}}-{{in.group2name}}.degs.txt")
write.table (out, degfile, quote=F, sep="\\t")

if ({{args.heatmap | R}}) {
	hmn  = min (nrow(out), {{args.hmn}})
	exps = expmatrix[row.names(out[1:hmn, ]), ]
	hmfile = file.path("{{degdir}}", "{{expfile | fn}}.{{in.group1name}}-{{in.group2name}}.heatmap.png")
	png (file = hmfile, res=300, width=2000, height=2000)
	heatmap(as.matrix(exps), margins = c({{args.hmmar}}), cexRow={{args.hmn}}/100)
	dev.off()
}

if ({{args.volplot | R}}) {
	dat = data.frame(out, n_log10_adj_pval = -c(log10(out$adj.P.Val)), col=ifelse(abs(out$logFC)>2 & out$adj.P.Val<0.01, 'A', ifelse(abs(out$logFC)>2, 'B', 'C')))
	a<-ggplot(dat, aes(x = logFC, y = n_log10_adj_pval, col=col))
	a<-a+ylab("-log10(adjusted P value)\n")
	a<-a+xlab("logFC")
	a<-a+theme_classic(base_size = 12)
	a<-a+theme(legend.position="none")
	a<-a+geom_point()
	a<-a+geom_vline(xintercept=c(2, -2))
	a<-a+geom_hline(yintercept = -log10(0.01))
	volfile = file.path("{{degdir}}", "{{expfile | fn}}.{{in.group1name}}-{{in.group2name}}.volcano.png")
	png (file = volfile, res=300, width=2000, height=2000)
	plot(a)
	dev.off()
}
"""

"""
@name:
	pRawCounts2
@description:
	Convert raw counts to another unit
@input:
	`expfile:file`: the expression matrix
	- rows are genes, columns are samples, if not use `args.transpose = True`
@output:
	`outfile:file`: the converted expression matrix
@args:
	`transpose`: transpose the input matrix? default: False
	`log2`:      whether to take log2? default: False
	`unit`:      convert to which unit? default: cpm (or rpkm, tmm)
	`header`:    whether input file has header? default: True
	`rownames`:  the index of the column as rownames. default: 1
	`glenfile`:  the gene length file, for RPKM
		- no head, row names are genes, have to be exact the same order and length as the rownames of expfile
@requires:
	[edgeR](https://bioconductor.org/packages/release/bioc/html/edger.html) if cpm or rpkm is chosen
	[coseq](https://rdrr.io/rforge/coseq/man/transform_RNAseq.html) if tmm is chosen
"""
pRawCounts2                = Proc(desc = 'Convert raw counts to another unit.')
pRawCounts2.input          = "expfile:file"
pRawCounts2.output         = "outfile:file:{{expfile | fn}}.{{args.unit}}.txt"
pRawCounts2.args.transpose = False
pRawCounts2.args.unit      = 'cpm'
pRawCounts2.args.header    = True
pRawCounts2.args.rownames  = 1
pRawCounts2.args.log2      = False
pRawCounts2.args.glenfile  = ''
pRawCounts2.lang           = "Rscript"
pRawCounts2.script         = """
data  = read.table ("{{expfile}}", sep="\\t", header={{args.header | R}}, row.names = {{args.rownames}}, check.names=F)
if ({{args.transpose | R}}) data = t (data)

if ("{{args.unit}}" == 'cpm') {
	library('edgeR')
	ret = cpm (data, log = {{args.log2 | R}})
} else if ("{{args.unit}}" == 'rpkm') {
	library('edgeR')
	genelen = read.table ("{{args.glenfile}}", header=F, row.names = 1, check.names = F)
	ret = rpkm (data, log = {{args.log2 | R}}, gene.length = as.vector(genelen))
} else {
	library('coseq')
	ret = transform_RNAseq(data, norm="TMM")
	ret = ret$normCounts
	if ({{args.log2 | R}})
		ret = log2(ret)
}

rnames = TRUE
cnames = TRUE
if ({{args.transpose | R}}) {
	rnames = {{args.header | R}}
	cnames = {{args.rownames}} == 1
} else {
	rnames = {{args.rownames}} == 1
	cnames = {{args.header | R}}
}

write.table (ret, "{{outfile}}", quote=F, row.names=rnames, col.names=cnames, sep="\\t")
"""