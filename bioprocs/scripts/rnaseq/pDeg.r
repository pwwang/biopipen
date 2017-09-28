

# get the exp data
ematrix    = read.table ("{{in.efile}}",  header=T, row.names = 1, check.names=F, sep="\t")
samples    = colnames(ematrix)

# get groups
gmatrix  = read.table ("{{in.gfile}}",  header=F, row.names = 1, check.names=F, sep="\t")
gfactor  = factor(gmatrix[,1,drop=T])
gfactor  = relevel(gfactor, gmatrix[1,1])
group1   = levels(gfactor)[1]
group2   = levels(gfactor)[2]
samples1 = rownames(gmatrix[gmatrix==group1,,drop=F])
samples2 = rownames(gmatrix[gmatrix==group2,,drop=F])
ematrix  = ematrix[, c(samples1, samples2), drop=F]

n1      = length(samples1)
n2      = length(samples2)
pval    = {{args.pval | lambda x: float(x)}}
filters = {% if args.filter | lambda x: isinstance(x, list) %}{{args.filter | Rvec}}{% else %}c({{args.filter}}){% endif %}
ematrix = ematrix[rowSums(ematrix>filters[1]) >= filters[2],,drop=F]
allfile = "{{out.outdir}}/{{in.efile | fn | fn}}-{{in.gfile | fn | fn}}.all.txt"

{% if args.paired %}
if (n1 != n2) {
	stop(paste("Flag paired is TRUE, but number of samples is different in two groups (", n1, n2, ").", sep=""))
}
{% endif %}

# detect DEGs using edgeR
{% if args.tool | lambda x: x == 'edger' %}
	# get model
	{% if args.paired %}
	pairs = c()
	group = c()
	for (i in 1:n1) {
		pairs = c(pairs, i, i)
		group = c(group, group1, group2)
	}
	pairs  = factor(pairs)
	group  = factor(group)
	design = model.matrix(~pairs + group)
	{% else %}
	group  = c(rep(group1, n1), rep(group2, n2))
	group  = factor(group)
	design = model.matrix(~group)
	{% endif %}

	library(edgeR)
	dge     = DGEList(counts = ematrix, group = group)
	dge     = dge[rowSums(cpm(dge)>filters[1]) >= filters[2], ]
	dge$samples$lib.size = colSums(dge$counts)
	dge     = calcNormFactors(dge)

	disp    = estimateDisp (dge, design)
	fit     = glmFit (disp, design)
	fit     = glmLRT (fit)

	allgene = topTags (fit, n=nrow(fit$table), p.value = 1)
	write.table (allgene$table, allfile, quote=F, sep="\t")

	degene  = allgene$table[allgene$table$PValue < {{args.pval}},,drop=F]
	write.table (degene, "{{out.outfile}}", quote=F, sep="\t")

	normedCounts = dge$counts
	alllogFC     = allgene$table$logFC
	allPval      = allgene$table$PValue


# detect DEGs using DESeq2
{% elif args.tool | lambda x: x == 'deseq2' %}
	library(DESeq2)
	# model
	{% if args.paired %}
	coldata           = data.frame(patients = factor(rep(1:n1, 2)), treats = factor(c(rep(group1, n1), rep(group2, n2))))
	coldata$treats    = relevel(coldata$treats, group1)
	rownames(coldata) = colnames(ematrix)
	dge               = DESeqDataSetFromMatrix(round(ematrix), coldata, design = ~patients + treats)
	{% else %}
	coldata           = data.frame(treats = factor(c(rep(group1, n1), rep(group2, n2))))
	coldata$treats    = relevel(coldata$treats, group1) # group1 frist
	rownames(coldata) = colnames(ematrix)
	dge               = DESeqDataSetFromMatrix(round(ematrix), coldata, design = ~treats)
	{% endif %}

	dge     = DESeq(dge)
	allgene = results(dge)
	write.table (allgene, allfile, quote=F, sep="\t")

	degene  = allgene[allgene$pvalue < {{args.pval}},,drop=F]
	write.table (degene, "{{out.outfile}}", quote=F, sep="\t")

	degene$logFC  = degene$log2FoldChange
	degene$PValue = degene$pvalue

	normedCounts  = log2(assay(dge) + 1)
	alllogFC      = allgene$log2FoldChange
	allPval       = allgene$pvalue

{% endif %}

# MDS plot 
{% if args.mdsplot %}
library(limma)
mdsplot = file.path ("{{out.outdir}}", "mdsplot.png")
do.call(png, c(list(filename=mdsplot), {{args.devpars | Rlist}}))
plotMDS(normedCounts)
dev.off()
{% endif %}

# Volcano
{% if args.volplot %}
{{ plotVolplot }}
volplot = file.path ("{{out.outdir}}", "volcano.png")
log2fc  = degene$logFC
fdr     = degene$PValue
plotVolplot(data.frame(logFC=log2fc, FDR=fdr), volplot, ggs = {{args.heatmapggs | Rlist}}, devpars = {{args.devpars | Rlist}})
{% endif %}

# heatmap
{% if args.heatmap %}
if (nrow(degene) < 2) {
	cat('pyppl.log.warning: Cannot generate heatmap as < 2 DEGs detected.', file = stderr())
} else {
	{{plotHeatmap}}
	hmap  = file.path ("{{out.outdir}}", "heatmap.png")
	ngene = nrow(degene)
	hmapn = {{args.heatmapn}}
	if (hmapn <= 1) hmapn = as.integer(hmapn * ngene)
	
	ngene = min(hmapn, nrow(degene))
	mat   = ematrix[1:ngene, ]
	plotHeatmap(mat, hmap, ggs = {{args.heatmapggs | Rlist}}, devpars = {{args.devpars | Rlist}})
}
{% endif %}

# maplots
{% if args.maplot %}
{{ plotMAplot }}
maplot = file.path("{{out.outdir}}", "maplot.png")
A = rowSums(ematrix) / ncol(ematrix)
M = alllogFC
threshold = allPval < {{args.pval}}

plotMAplot(data.frame(A, M, threshold), maplot, ggs = {{args.maplotggs | Rlist}}, devpars = list(height=2000,res=300,width=2000))
		
{% endif %}