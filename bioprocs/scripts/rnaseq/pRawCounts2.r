{{rimport}}('plot.r')

infile  = {{in.infile | R}}
outfile = {{out.outfile | R}}
outdir  = {{out.outdir | R}}
prefix  = {{in.infile | fn2 | R}}
unit    = {{args.unit | R}}
refgene = {{args.refgene | R}}
hmrows  = {{args.hmrows | R}}
plot    = {{args.plot | R}}
ggs     = {{args.ggs | R}}
devpars = {{args.devpars | R}}

data    = read.table.nodup (infile, sep="\t", header=T, row.names=1, check.names=F)
rnames  = rownames(data)

if (unit == 'cpm') {
	library('edgeR')
	ret = cpm(data, log = F)
} else if (unit == 'fpkm' || unit == 'rpkm') {
	library('edgeR')
	genelen = read.table (refgene, header=F, row.names = NULL, check.names = F, sep = "\t")
	genelen = apply(genelen, 1, function(row) {
		gene = unlist(strsplit(row[9], ' ', fixed=T))[2]
		gene = substr(gene, 1, nchar(gene) - 1)
		len  = as.numeric(row[5]) - as.numeric(row[4])
		c(gene, len)
	})
	genelen = t(genelen)
	genelen = as.vector(genelen[which(genelen[,1] %in% rnames),2,drop=T])
	ret = rpkm (data, log = F, gene.length = as.numeric(genelen))
} else if (unit == 'tmm') {
	library('coseq')
	ret = transform_RNAseq(data, norm="TMM")
	ret = ret$normCounts
} else {
	stop(paste('Unknown expression unit:', unit))
}

write.table (round(ret, 3), outfile, quote=F, row.names=T, col.names=T, sep="\t")

exp = log2(ret + 1)

if (plot$boxplot) {
	bpfile = file.path(outdir, paste0(prefix, '.boxplot.png'))
	plot.boxplot(exp, bpfile, stack = T, devpars = devpars, ggs = ggs$boxplot)
}

if (plot$heatmap) {
	hmfile = file.path(outdir, paste0(prefix, ".heatmap.png"))
	hmexp  = if (nrow(exp) > hmrows) exp[sample(nrow(exp),size=hmrows),] else exp
	plot.heatmap(hmexp, hmfile, devpars = devpars, ggs = ggs$heatmap)
}

if (plot$histogram) {
	histfile = file.path(outdir, paste0(prefix, ".histo.png"))
	plot.histo(stack(as.data.frame(exp)), histfile, devpars = devpars, ggs = ggs$histogram)
}
