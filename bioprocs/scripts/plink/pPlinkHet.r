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
	het   = T,
	out   = output
)
cmd = sprintf("%s %s 1>&2", plink, cmdargs(params, equal = ' '))
runcmd(cmd)

phet = read.table(paste0(output, '.het'), header = T, row.names = NULL, check.names = F)
het = data.frame(Het = 1 - phet[, "O(HOM)"]/phet[, "N(NM)"])
rownames(het) = paste(phet$FID, phet$IID, sep = "\t")
het.mean = mean(het$Het, na.rm = T)
het.sd   = sd(het$Het, na.rm = T)
het.fail = rownames(het[!is.na(het$Het) & (het$Het < het.mean-cutoff*het.sd | het$Het > het.mean+cutoff*het.sd),, drop = F])
writeLines(het.fail, con = file(paste0(output, '.het.fail')))

{% if args.plot %}
het$Status = "Pass"
het[het.fail, "Status"] = "Fail"

plot.histo(
	data     = het,
	plotfile = paste0(output, '.het.png'),
	params   = list(aes(fill=Status), bins = 50),
	ggs      = list(
		xlab       = list("Sample Heterozygosity"),
		ylab       = list("Count"),
		geom_vline = list(xintercept = c(het.mean-cutoff*het.sd, het.mean+cutoff*het.sd), color = "red", linetype="dashed"),
		geom_vline = list(xintercept = het.mean, color = "blue", linetype="dashed"),
		theme      = list(legend.position = "none"),
		geom_text  = list(
			aes(x = het.mean-cutoff*het.sd, y = Inf, label = sprintf('mean - %ssd (%.3f)', cutoff, het.mean - cutoff*het.sd)),
			colour="red", angle=90, vjust = 1.2, hjust = 1.2
		),
		geom_text  = list(
			aes(x = het.mean+cutoff*het.sd, y = Inf, label = sprintf('mean + %ssd (%.3f)', cutoff, het.mean + cutoff*het.sd)),
			colour="red", angle=90, vjust = 1.2, hjust = 1.2
		),
		geom_text  = list(
			aes(x = het.mean, y = Inf, label = sprintf('mean (%.3f)', het.mean)),
			colour="blue", vjust = 1.5, hjust = -.1
		)
	),
	devpars = devpars
)
{% endif %}