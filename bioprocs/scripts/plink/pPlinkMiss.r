library(methods)
{{rimport}}('__init__.r', 'plot.r')

indir  = {{i.indir | R}}
outdir = {{o.outdir | R}}

plink    = {{args.plink | R}}
samplecr = {{args.samplecr | R}}
snpcr    = {{args.snpcr | R}}
devpars  = {{args.devpars | R}}

bedfile = Sys.glob(file.path(indir, '*.bed'))
input   = tools::file_path_sans_ext(bedfile)
output  = file.path(outdir, basename(input))

params = list(
	bfile   = input,
	missing = T,
	out     = output
)
cmd = sprintf("%s %s 1>&2", plink, cmdargs(params, equal = ' '))
runcmd(cmd)

imiss = read.table(paste0(output, '.imiss'), header = T, row.names = NULL, check.names = F)
callrate.sample = data.frame(Callrate = 1-imiss$F_MISS)
rownames(callrate.sample) = paste(imiss$FID, imiss$IID, sep = "\t")
callrate.sample.fail = rownames(callrate.sample[callrate.sample$Callrate < samplecr, , drop = F])
writeLines(callrate.sample.fail, con = file(paste0(output, '.samplecr.fail')))

{% if args.plot %}
callrate.sample$Status = "Pass"
callrate.sample[callrate.sample.fail, "Status"] = "Fail"
plot.histo(
	data     = callrate.sample,
	plotfile = paste0(output, '.samplecr.png'),
	params   = list(aes(fill=Status), bins = 50),
	ggs      = list(
		xlab       = list("Sample Call Rate"),
		ylab       = list("Count"),
		geom_vline = list(xintercept = samplecr, color = "red", linetype="dashed"),
		theme      = list(legend.position = "none"),
		geom_text  = list(
			aes(x = samplecr, y = Inf, label = samplecr),
			colour="red", angle=90, vjust = 1.2, hjust = 1.2
		)
	)
)
{% endif %}

lmiss = read.table(paste0(output, '.lmiss'), header = T, row.names = NULL, check.names = F)
lmiss$Callrate = 1-lmiss$F_MISS
callrate.snp.fail = lmiss[which(lmiss$Callrate < snpcr), 'SNP', drop = F]
write.table(callrate.snp.fail, paste0(output, '.snpcr.fail'), row.names = F, col.names = F, sep = "\t", quote = F)

{% if args.plot %}
lmiss$Status = "Pass"
lmiss[which(lmiss$Callrate < snpcr), "Status"] = "Fail"
plot.histo(
	data     = lmiss,
	plotfile = paste0(output, '.snpcr.png'),
	x        = 'Callrate',
	params   = list(aes(fill=Status), bins = 50),
	ggs      = list(
		xlab       = list("SNP Call Rate"),
		ylab       = list("Count"),
		geom_vline = list(xintercept = snpcr, color = "red", linetype="dashed"),
		theme      = list(legend.position = "none"),
		geom_text  = list(
			aes(x = snpcr, y = Inf, label = snpcr),
			colour="red", angle=90, vjust = 1.2, hjust = 1.2
		)
	),
	devpars = devpars
)
{% endif %}


