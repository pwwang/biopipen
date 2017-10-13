{{mem2}}
{{params2CmdArgs}}

params = {{args.params | Rlist}}
params$i = {{in.infile | quote}}
params$o = {{out.outfile | quote}}

cmd = paste('{{args.bamstats}}', params2CmdArgs(params), sep = ' ')

{% if args.plot %}
{{plotHist}}
{{plotBoxplot}}
{{pollingFirst}}
pollingFirst ({{proc.workdir | quote}}, {{proc.size}}, {{job.index}}, cmd, "bamstats.done")

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

write ("Plotting average coverages ...", stderr())
write.table (means, "{{out.outdir}}/avgCoverage.txt", quote=F, sep="\\t")
plotHist (means, "{{out.outdir}}/avgCoverage.png", ggs = {{args.histplotggs | Rlist}}, devpars = {{args.devpars | Rlist}})

# plot chromosomes
write ("Plotting chromosome coverages ...", stderr())
png ("{{out.outdir}}/chrCoverage.png")
plotBoxplot(t(chrs), "{{out.outdir}}/chrCoverage.png", ggs = {{args.boxplotggs | Rlist}}, devpars = {{args.devpars | Rlist}})
dev.off()

{% endif %}

{% else %}
{{runcmd}}
runcmd (cmd)
{% endif %}

