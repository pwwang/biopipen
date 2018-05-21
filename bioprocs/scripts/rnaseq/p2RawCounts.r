{{rimport}}('plot.r')

infile  = {{in.infile | R}}
outfile = {{out.outfile | R}}
outdir  = {{out.outdir | R}}
prefix  = {{in.infile | fn2 | R}}
unit    = {{args.unit | R}}
refgene = {{args.refgene | R}}
nreads  = {{args.nreads | R}}
hmrows  = {{args.hmrows | R}}
plot    = {{args.plot | R}}
ggs     = {{args.ggs | R}}
devpars = {{args.devpars | R}}

data    = read.table.nodup (infile, sep="\t", header=T, row.names=1, check.names=F)
rnames  = rownames(data)

getGenelen = function() {
	gdata  = read.table(refgene, sep="\t", header=F, row.names=NULL, check.names=F)
	gnames = strsplit(as.vector(gdata[, ncol(gdata)]), ";", fixed=T)
	gnames = matrix(unlist(gnames), nrow=length(gnames[[1]]))[1,]
	gnames = strsplit(gnames, " ", fixed=T)
	gnames = matrix(unlist(gnames), nrow=length(gnames[[1]]))
	gnames = gnames[nrow(gnames),]
	gnames = make.unique(gnames)
	glen   = matrix(gdata[,5] - gdata[,4], ncol = 1)/1000
	colnames(glen) = c("GLEN")
	rownames(glen) = gnames
	return (glen)
}

if (unit == 'fpkm' || unit == 'rpkm') {
	# convert to total 50,000,000 reads per sample
	glen = getGenelen()
	rn   = intersect(rownames(glen), rnames)
	ret  = data[rn,,drop=F] * glen[rn,,drop=F]
	ssum = colSums(ret)
	ret  = nreads * ret / ssum
} else if (unit == 'tpm') {
	glen = getGenelen()
	rn   = intersect(rownames(glen), rnames)
	ret  = data[rn,,drop=F]
	ssum = colSums(ret)
	ret  = nreads * ret / ssum
	ret  = ret * glen[rn,,drop=F] / sum(glen[rn,,drop=T])
}

write.table (round(ret), outfile, quote=F, row.names=T, col.names=T, sep="\t")

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

