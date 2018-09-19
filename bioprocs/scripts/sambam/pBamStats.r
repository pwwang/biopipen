{{rimport}}('__init__.r', 'poll.r')

params = {{args.params | Rlist}}
params$i = {{i.infile | quote}}
params$o = {{o.outfile | quote}}
{% if args.feature != 'wgs' %}
params$f = {{args.feature | quote}}
{% endif %}

cmd = paste('{{args.bamstats}}', mem2({{args.mem | quote}}, 'java'), cmdargs(params), sep = ' ')
poll = Poll({{proc.workdir | quote}}, {{proc.size}}, {{job.index}})

{% if args.plot %}
{{rimport}}('plot.r')

{% if job.index == 0 %}
runcmd(cmd)
{% endif %}
poll$non1st(cmd)

{% if job.index == 0 %}
##### start plotting

bsfiles = Sys.glob("{{proc.workdir}}/*/output/*/*.stat.txt")
means   = matrix(ncol=1, nrow=length(bsfiles))
chrs    = NULL

rnames  = make.unique(unlist(lapply(bsfiles, function(x){ x=basename(x); return (substr(x, 1, nchar(x)-9)) })), sep='_')
rownames(means) = rnames
colnames(means) = c("Average_coverage")
for (i in 1:length(bsfiles)) {
	bsfile = bsfiles[i]
	if (file.info(bsfile)$size == 0) {
		logger ('No content in file:', bsfile)
		next
	}
	logger ("Reading", bsfile, "...")
	sample = rnames[i]
	stat   = read.table (bsfile, sep="", header=T, check.names=F, row.names=1)
	#stat   = stat[chrs2, ]
	#stat[, "N"] = as.numeric(gsub(",", "", stat[, "N"]))
	m      = stat[which(stat[, "mean"] > {{args.cutoff}} & stat[, "mean"] < {{args.cap}}), "mean", drop=F]
	#N      = stat[rownames(m), "N", drop = F]
	#means[sample, 1] = sum(N * m)/sum(N)
	means[sample, 1] = mean(m[, "mean", drop=T])
	col2in = m[order(m[, "mean"], decreasing = T), "mean", drop=F][1:(2*{{args.nfeats}}), , drop=F]
	colnames(col2in) = sample
	if (is.null(chrs)) {
		chrs = col2in
	} else {
		chrs = cbindfill(chrs, col2in)
	}
}
chrs = chrs[order(rowMeans(chrs), decreasing = T),,drop=F][1:min({{args.nfeats}}, nrow(chrs)),,drop=F]

logger("Plotting average coverages ...")
write.table (means, "{{out.outdir}}/avgCoverage.txt", quote=F, sep="\t")
plot.histo (means, "{{out.outdir}}/avgCoverage.png", ggs = {{args.histplotggs | Rlist}}, devpars = {{args.devpars | Rlist}})

# plot chromosomes
logger ("Plotting feature coverages ...")
plot.boxplot(t(chrs), "{{out.outdir}}/featureCoverage.png", ggs = {{args.boxplotggs | Rlist}}, devpars = {{args.devpars | Rlist}})

{% endif %} # end if job.index

{% else %}
runcmd (cmd)
{% endif %} # end if args.plot
