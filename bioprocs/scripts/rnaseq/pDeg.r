

# get the exp data
ematrix   = read.table ("{{in.efile}}",  header=T, row.names = 1, check.names=F, sep="\t")

# get groups
{{ txtSampleinfo }}
sampleinfo = txtSampleinfo("{{in.gfile}}")

gfactor   = factor(sampleinfo$Group)
gfactor   = relevel(gfactor, as.vector(sampleinfo$Group[1]))
group1    = levels(gfactor)[1]
group2    = levels(gfactor)[2]
samples1  = rownames(sampleinfo[which(sampleinfo$Group == group1), , drop=F])
samples2  = rownames(sampleinfo[which(sampleinfo$Group == group2), , drop=F])
samples   = c(samples1, samples2)
ematrix   = ematrix[, samples, drop=F]

n1      = length(samples1)
n2      = length(samples2)
pval    = {{args.pval | lambda x: float(x)}}
filters = {% if args.filter | lambda x: isinstance(x, list) %}{{args.filter | Rvec}}{% else %}c({{args.filter}}){% endif %}
ematrix = ematrix[rowSums(ematrix>filters[1]) >= filters[2],,drop=F]
allfile = "{{out.outdir}}/{{in.efile | fn | fn}}-{{in.gfile | fn | fn}}.all.txt"

if ("Patient" %in% colnames(sampleinfo) && n1 != n2) {
	stop(paste("Paired samples indicated, but number of samples is different in two groups (", n1, n2, ").", sep=""))
}

# detect DEGs using edgeR
{% if args.tool | lambda x: x == 'edger' %}
	# get model
	if ("Patient" %in% colnames(sampleinfo)) {
		pairs  = factor(sampleinfo[samples, "Patient"])
		group  = factor(sampleinfo[samples, "Group"])
		group  = relevel(group, group2)
		design = model.matrix(~pairs + group)
	} else {
		group  = factor(sampleinfo[samples, "Group"])
		group  = relevel(group, group2)
		design = model.matrix(~group)
	}

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

	degene  = allgene$table[allgene$table$PValue < pval,,drop=F]
	write.table (degene, "{{out.outfile}}", quote=F, sep="\t")

	normedCounts = dge$counts
	alllogFC     = allgene$table$logFC
	allPval      = allgene$table$PValue


# detect DEGs using DESeq2
{% elif args.tool | lambda x: x == 'deseq2' %}
	library(DESeq2)
	# model
	if ("Patient" %in% colnames(sampleinfo)) {
		coldata           = data.frame(patients = factor(sampleinfo[samples, 'Patient']), treats = factor(sampleinfo[samples, 'Group']))
		coldata$treats    = relevel(coldata$treats, group2)
		rownames(coldata) = samples
		dge               = DESeqDataSetFromMatrix(round(ematrix), coldata, design = ~patients + treats)
	} else {
		coldata           = data.frame(treats = factor(sampleinfo[samples, 'Group']))
		coldata$treats    = relevel(coldata$treats, group2) 
		rownames(coldata) = samples
		dge               = DESeqDataSetFromMatrix(round(ematrix), coldata, design = ~treats)
	}

	dge     = DESeq(dge)
	allgene = results(dge)
	write.table (allgene, allfile, quote=F, sep="\t")

	degene  = allgene[allgene$pvalue < pval,,drop=F]
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
plotMDS(normedCounts, col=c(rep("red", n1), rep("blue", n2)))
dev.off()
{% endif %}

# Volcano
{% if args.volplot %}
{{ plotVolplot }}
volplot = file.path ("{{out.outdir}}", "volcano.png")
log2fc  = alllogFC
fdr     = allPval
fdrcut  = rep(pval, length(allPval))
plotVolplot(data.frame(logFC=log2fc, FDR=fdr, FDRCut=fdrcut), volplot, ggs = {{args.heatmapggs | Rlist}}, devpars = {{args.devpars | Rlist}})
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
	mat   = ematrix[rownames(degene[1:ngene,,drop=F]), ]
	plotHeatmap(mat, hmap, ggs = {{args.heatmapggs | Rlist}}, devpars = {{args.devpars | Rlist}})
}
{% endif %}

# maplots
{% if args.maplot %}
{{ plotMAplot }}
maplot = file.path("{{out.outdir}}", "maplot.png")
A = rowSums(normedCounts) / ncol(normedCounts)
M = alllogFC
threshold = allPval < pval

plotMAplot(data.frame(A, M, threshold), maplot, ggs = {{args.maplotggs | Rlist}}, devpars = list(height=2000,res=300,width=2000))
		
{% endif %}