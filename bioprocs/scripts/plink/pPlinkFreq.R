library(methods)
{{'__init__.R', 'plot.R' | rimport}}

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
	freq = T,
	out   = output
)
shell$load_config(plink = plink)
shell$plink(params)$fg
#//cmd = sprintf("%s %s 1>&2", plink, cmdargs(params, equal = ' '))
#//runcmd(cmd)

freq = read.table(paste0(output, '.frq'), header = T, row.names = NULL, check.names = F)
freq.fail = freq[which(freq$MAF < cutoff), 'SNP', drop = F]

write.table(freq.fail, paste0(output, '.frq.fail'),
            col.names = F, row.names = F, sep = "\t", quote = F)

{% if args.plot %}
freq$Status = "Pass"
freq[which(freq$SNP %in% freq.fail$SNP), "Status"] = "Fail"

plot.histo(
	data     = freq,
	x        = 'MAF',
	plotfile = paste0(output, '.frq.png'),
	params   = list(aes(fill=Status)),
	ggs      = list(
		xlab       = list("MAF"),
		ylab       = list("Count"),
		geom_vline = list(xintercept = cutoff, color = "red", linetype="dashed"),
		theme      = list(legend.position = "none"),
		geom_text  = list(
			aes(x = cutoff, y = Inf, label = cutoff),
			colour="red", angle=90, vjust = 1.2, hjust = 1.2
		)
	),
	devpars = devpars
)
{% endif %}
