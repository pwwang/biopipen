from pyppl import proc

"""
@name:
	pCallByLimmaFromMatrix
@description:
	Call DEG from expressoin matrix, where column names must in accordant order of <group>
@input:
	`matfile:file`: the expression matrix
	`group1`:       columns of group1 (separated by comma)
	`group2`:       columns of group2 (separated by comma)
	`group1name`:   the name of group1
	`group2name`:   the name of group2
@output:
	`degfile:file`: the output file containing DEGs
@args:
	`pval`: the cutoff of DEGs (default: .05)
@requires:
	[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
"""
pCallByLimmaFromMatrix = proc ()
pCallByLimmaFromMatrix.input     = "matfile:file, group1, group2, group1name, group2name"
pCallByLimmaFromMatrix.output    = "degfile:file:{{matfile.fn}}.deg"
pCallByLimmaFromMatrix.args      = {'pval': 0.05}
pCallByLimmaFromMatrix.defaultSh = "Rscript"
pCallByLimmaFromMatrix.script    = """
library('limma')
expdata = read.table ("{{matfile}}", sep="\\t", header=T, row.names=1)
group1  = strsplit("{{group1}}", ",")[[1]]
group2  = strsplit("{{group2}}", ",")[[1]]
g1data  = expdata[, which(names(expdata) %in% group1)]
g2data  = expdata[, which(names(expdata) %in% group2)]

group   = c(rep("{{group1name}}", length(group1)), rep("{{group2name}}", length(group2)))
design  = model.matrix (~group)
fit     = lmFit (cbind(g1data, g2data), design)
fit     = eBayes(fit)
ret     = topTable(fit, n=length(fit$p.value))
out     = ret [ret$P.Value < {{proc.args.pval}}, ]
deg     = paste (rownames(out), collapse=',')
opv     = paste (format(out$P.Value, digits=3, scientific=T), collapse=',')
write (paste(deg, opv, sep="\\t"), file = "{{degfile}}", append=T)
"""

"""
@name:
	pCallByLimmaFromFiles
@description:
	Call DEG from expression files
@input:
	`expdir:file`:  the directory containing expression files
	`group1`:       columns of group1 (separated by comma)
	`group2`:       columns of group2 (separated by comma)
	`group1name`:   the name of group1
	`group2name`:   the name of group2   
@output:
	`degfile:file`: the output file containing DEGs
@args:
	`pval`: the cutoff of DEGs (default: .05)
	`paired`: whether the samples are paired
@requires:
	[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
"""
pCallByLimmaFromFiles = proc ()
pCallByLimmaFromFiles.input     = "expdir:file, group1, group2, group1name, group2name"
pCallByLimmaFromFiles.output    = "degfile:file:{{expdir.fn}}.deg"
pCallByLimmaFromFiles.args      = {'pval': 0.05, 'paired': False}
pCallByLimmaFromFiles.defaultSh = "Rscript"
pCallByLimmaFromFiles.script    = """
library('limma')
group1  = strsplit("{{group1}}", ",")[[1]]
group2  = strsplit("{{group2}}", ",")[[1]]

expmatrix = matrix()
for (i in 1:length (group1)) {
	file = paste ("{{expdir}}/", group1[i], sep="")
	if (grepl ('.gz$', group1[i])) {
		file = gzfile(file)
	}
	tmp  = read.table (file, sep="\\t", header=F, row.names = 1)
	if (i==1) {
		expmatrix = as.matrix(tmp)
	}
	else {
		expmatrix = cbind (expmatrix, tmp)
	}
}

for (i in 1:length (group2)) {
	file = paste ("{{expdir}}/", group2[i], sep="")
	if (grepl ('.gz$', group2[i])) {
		file = gzfile(file)
	}
	tmp  = read.table (file, sep="\\t", header=F, row.names = 1)
	expmatrix = cbind (expmatrix, tmp)
}
pairs = vector(mode="numeric")
group = vector(mode="character")
ng1   = length(group1)
ng2   = length(group2)
if ({{proc.args.paired | str(_).upper()}}) {
	for (i in 1:ng1) {
		pairs = c(pairs, i, i)
		group = c(group, "{{group1name}}", "{{group2name}}")
	}
	pairs = factor(pairs)
	group = factor(group)
	design = model.matrix(~pairs+group)
} else {
	group   = c(rep("{{group1name}}", ng1), rep("{{group2name}}", ng2))
	group = factor(group)
	design = model.matrix(~group)
}
fit    = lmFit (expmatrix, design)
fit    = eBayes(fit)
if ({{proc.args.paired | str(_).upper()}}) {
	ret    = topTable(fit, n=length(fit$p.value), coef=paste("group", "{{group1name}}", sep=""))
} else {
	ret    = topTable(fit, n=length(fit$p.value))
}
out    = ret [ret$P.Value < {{proc.args.pval}}, ]
write.table (out, "{{degfile}}", quote=F, sep="\\t")
"""

"""
@name:
	pCallByEdgeRFromFiles
@description:
	Call DEG from expression files
@input:
	`expdir:file`:  the directory containing expression files
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
	`paired`:  whether the samples are paired
	`bcvplot`: whether to plot biological coefficient of variation, default: True
	`displot`: whether to plot biological coefficient of variation, default: True
	`fcplot`:  whether to plot fold changes, default: True
@requires:
	[edgeR](https://bioconductor.org/packages/release/bioc/html/edger.html)
"""
pCallByEdgeRFromFiles = proc ()
pCallByEdgeRFromFiles.input     = "expdir:file, group1, group2, group1name, group2name"
pCallByEdgeRFromFiles.output    = "degdir:dir:{{expdir.fn}}.deg"
pCallByEdgeRFromFiles.args      = {'filter': "1,2", 'pval': 0.05, 'paired': False, 'bcvplot': True, 'displot': True, 'fcplot': True}
pCallByEdgeRFromFiles.defaultSh = "Rscript"
pCallByEdgeRFromFiles.script    = """
library('edgeR')
group1  = strsplit("{{group1}}", ",")[[1]]
group2  = strsplit("{{group2}}", ",")[[1]]

expmatrix = matrix()
for (i in 1:length (group1)) {
	file = paste ("{{expdir}}/", group1[i], sep="")
	if (grepl ('.gz$', group1[i])) {
		file = gzfile(file)
	}
	tmp  = read.table (file, sep="\\t", header=F, row.names = 1)
	colnames (tmp) = c (group1[i])
	if (i==1) {
		expmatrix = as.matrix(tmp)
	}
	else {
		expmatrix = cbind (expmatrix, tmp)
	}
}

for (i in 1:length (group2)) {
	file = paste ("{{expdir}}/", group2[i], sep="")
	if (grepl ('.gz$', group2[i])) {
		file = gzfile(file)
	}
	tmp  = read.table (file, sep="\\t", header=F, row.names = 1)
	colnames (tmp) = c (group2[i])
	expmatrix = cbind (expmatrix, tmp)
}
pairs = vector(mode="numeric")
group = vector(mode="character")
ng1   = length(group1)
ng2   = length(group2)
if ({{proc.args.paired | str(_).upper()}}) {
	for (i in 1:ng1) {
		pairs = c(pairs, i, i)
		group = c(group, "{{group1name}}", "{{group2name}}")
	}
	pairs  = factor(pairs)
	group  = factor(group)
	design = model.matrix(~pairs+group)
} else {
	group  = c(rep("{{group1name}}", ng1), rep("{{group2name}}", ng2))
	group  = factor(group)
	design = model.matrix(~group)
}

# filter
dobj   = DGEList(counts=expmatrix, group=group)
filter = noquote(unlist(strsplit("{{proc.args.filter}}", ",")))
fX     = as.numeric (filter[1])
fY     = as.numeric (filter[2])
dobj   = dobj[rowSums(cpm(dobj)>fX) >= fY, ]
dobj$samples$lib.size = colSums(dobj$counts)

# normalize
dobj = calcNormFactors(dobj, method="TMM")

if ({{proc.args.bcvplot | str(_).upper()}}) {
	bcvplot = file.path ("{{degdir}}", "bcvplot.png")
	png (file=bcvplot)
	plotMDS (dobj, method="bcv", col=as.numeric(dobj$samples$group))
	legend("bottomleft", as.character(unique(dobj$samples$group)), col=2:1, pch=20)
	dev.off()
}

disp <- estimateDisp (dobj, design)
if ({{proc.args.displot | str(_).upper()}}) {
	displot = file.path ("{{degdir}}", "displot.png")
	png (file=displot)
	plotBCV (disp)
	dev.off()
}

fit    = glmFit (disp, design)
fit    = glmLRT (fit)

if ({{proc.args.fcplot | str(_).upper()}}) {
	deg    = decideTestsDGE(fit, p.value = {{proc.args.pval}})
	fcplot = file.path ("{{degdir}}", "fcplot.png")
	png (file=fcplot)
	tags = rownames(disp)[as.logical(deg)]
	plotSmear (fit, de.tags=tags)
	abline(h = c(-1, 1), col = "blue")
	dev.off()
}

out    = topTags (fit, n=nrow(fit$table), p.value = {{proc.args.pval}})
write.table (out$table, file.path("{{degdir}}", "degs.txt"), quote=F, sep="\\t")
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
	`unit`:      convert to which unit? default: cpm (or rpkm)
	`header`:    whether input file has header? default: True
	`rownames`:  the index of the column as rownames. default: 1
	`glenfile`:  the gene length file, for RPKM
	- no head, row names are genes, have to be exact the same order and length as the rownames of expfile
@requires:
	[edgeR](https://bioconductor.org/packages/release/bioc/html/edger.html)
"""
pRawCounts2 = proc ()
pRawCounts2.input     = "expfile:file"
pRawCounts2.output    = "outfile:file:{{expfile.fn}}.{{proc.args.unit}}.txt"
pRawCounts2.args      = {'transpose': False, 'unit': 'cpm', 'header': True, 'rownames': 1, 'log2': False, 'glenfile': ''}
pRawCounts2.defaultSh = "Rscript"
pRawCounts2.script    = """
library('edgeR')

data  = read.table ("{{expfile}}", sep="\\t", header={{proc.args.header | str(_).upper()}}, row.names = {{proc.args.rownames}}, check.names=F)
if ({{proc.args.transpose | str(_).upper()}}) data = t (data)

if ("{{proc.args.unit}}" == 'cpm') {
	ret = cpm (data, log = {{proc.args.log2 | str(_).upper()}})
} else {
	genelen = read.table ("{{proc.args.glenfile}}", header=F, row.names = 1, check.names = F)
	ret = rpkm (data, log = {{proc.args.log2 | str(_).upper()}}, gene.length = as.vector(genelen))
}

rnames = TRUE
cnames = TRUE
if ({{proc.args.transpose | str(_).upper()}}) {
	rnames = {{proc.args.header | str(_).upper()}}
	cnames = {{proc.args.rownames}} == 1
} else {
	rnames = {{proc.args.rownames}} == 1
	cnames = {{proc.args.header | str(_).upper()}}
}

write.table (ret, "{{outfile}}", quote=F, row.names=rnames, col.names=cnames, sep="\\t")
"""