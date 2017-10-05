

# get the exp data
ematrix    = read.table ("{{in.efile}}",  header=T, row.names = 1, check.names=F, sep="\t")

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

# detect DEGs using limma
{% if args.tool | lambda x: x == 'limma' %}
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

	library(limma)
	fit     = lmFit(ematrix, design)
	fit     = eBayes(fit)

	allgene = topTable(fit, n=nrow(ematrix), adjust="BH", coef = paste("group", group1, sep=""))
	write.table (allgene, allfile, quote=F, sep="\t")

	degene  = allgene[allgene$P.Value < pval,,drop=F]
	write.table (degene, "{{out.outfile}}", quote=F, sep="\t")

	normedCounts = ematrix
	alllogFC     = allgene$logFC
	allPval      = allgene$P.Value

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
A = rowSums(ematrix) / ncol(ematrix)
M = alllogFC
threshold = allPval < pval

plotMAplot(data.frame(A, M, threshold), maplot, ggs = {{args.maplotggs | Rlist}}, devpars = list(height=2000,res=300,width=2000))
		
{% endif %}