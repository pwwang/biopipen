library(methods)
{{rimport}}('__init__.r', 'plot.r')

indir  = {{i.indir | R}}
outdir = {{o.outdir | R}}

plink   = {{args.plink | R}}
cutoff  = {{args.cutoff | R}}
devpars = {{args.devpars | R}}

bedfile = Sys.glob(file.path(indir, '*.bed'))
input   = tools::file_path_sans_ext(bedfile)
output  = file.path(outdir, basename(input))

params = list(
	bfile = input,
	hardy = T,
	out   = output
)
cmd = sprintf("%s %s 1>&2", plink, cmdargs(params, equal = ' '))
runcmd(cmd)

hardy = read.table(paste0(output, '.hwe'), header = T, row.names = NULL, check.names = F)
hardy.fail = hardy[which(hardy$P < cutoff), 'SNP', drop = F]
write.table(hardy.fail, paste0(output, '.hardy.fail'), col.names = F, row.names = F, sep = "\t", quote = F)

{% if args.plot %}
hardy$Pval   = -log10(hardy$P)
hardy$Status = "Pass"
hardy[which(hardy$SNP %in% hardy.fail$SNP), "Status"] = "Fail"

plot.histo(
	data     = hardy,
	x        = 'Pval',
	plotfile = paste0(output, '.hardy.png'),
	params   = list(aes(fill=Status), bins = 50),
	ggs      = list(
		xlab       = list("-log10(HWE p-value)"),
		ylab       = list("Count"),
		geom_vline = list(xintercept = -log10(cutoff), color = "red", linetype="dashed"),
		theme      = list(legend.position = "none"),
		geom_text  = list(
			aes(x = -log10(cutoff), y = Inf, label = cutoff),
			colour="red", angle=90, vjust = 1.2, hjust = 1.2
		)
	),
	devpars = devpars
)
{% endif %}