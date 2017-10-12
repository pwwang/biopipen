{{mem2}}

cmd = '{{args.bamstats}} -i "{{in.infile}}" -o "{{out.outfile}}" {{args.params}}'

{% if args.plot %}
{{pollingAll}}
pollingAll ({{proc.workdir | quote}}, {{proc.size}}, {{job.index}}, cmd, "bamstats.done")

{% if job.index | lambda x: x == 0 %}
##### start plotting

bsfiles = Sys.glob("{{proc.workdir}}/*/output/*/*.stat.txt")
means   = matrix(ncol=1, nrow=length(bsfiles))
chrs    = NULL

rnames  = make.unique(unlist(lapply(bsfiles, function(x){ x=basename(x); return (substr(x, 1, nchar(x)-9)) })), sep='_')
rownames(means) = rnames
colnames(means) = c("Average coverage")
for (i in 1:length(bsfiles)) {
	bsfile = bsfiles[i]
	write (paste("Reading", bsfile, "...", sep=" "), stderr())
	sample = rnames[i]
	stat   = read.table (bsfile, sep="", header=T, check.names=F, row.names=1)
	#stat   = stat[chrs2, ]
	stat[, "N"] = as.numeric(gsub(",", "", stat[, "N"]))
	means[sample, 1] = sum(stat[, "N"] * stat[, "mean"])/sum(stat[, "N"])
	col2in = stat[, "mean", drop=F]
	colnames(col2in) = sample
	if (is.null(chrs)) {
		chrs = col2in
	} else {
		chrs = cbind(chrs, col2in)
	}
}
# plot average coverage
plotFreq = function (obj, figure, xlab, ylab="Frequency") {
	png (file=figure)
	h = hist (obj, freq=T, xlab=xlab, ylab=ylab, col="gray", main=paste(xlab, "distribution", sep=" "), axes=F)
	minb = min(h$breaks)
	maxb = max(h$breaks)
	maxc = max(h$counts)
	lenb = length(h$breaks)
	stpb = (maxb-minb)/(lenb-1)
	axis(1, pos=0, labels=T, at=seq(minb,maxb,stpb))
	lab0 = floor(log10(maxc))
	stpc = ceiling(maxc/(10**lab0)) * (10 ** (lab0-1))
	axis(2, pos=minb, labels=T, at=seq(0, maxc, stpc))
	dev.off()
}
write ("Plotting average coverages ...", stderr())
write.table (means, "{{out.outdir}}/avgCoverage.txt", quote=F, sep="\\t")
plotFreq (means, "{{out.outdir}}/avgCoverage.png", xlab="Average coverage")

# plot chromosomes
write ("Plotting chromosome coverages ...", stderr())
png ("{{out.outdir}}/chrCoverage.png")
#colnames(chrs) = NULL # just show index
boxplot(t(chrs), ylab="Coverage", las=2)
dev.off()

{% endif %}

{% else %}
{{runcmd}}
runcmd (cmd)
{% endif %}

