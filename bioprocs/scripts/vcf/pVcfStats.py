
from bioprocs.utils import shell2 as shell

infile   = {{i.infile | quote}}
config   = {{i.config | quote}}
outdir   = {{o.outdir | quote}}
Rscript  = {{args.Rscript | quote}}
formula  = {{args.formula | repr}}
title    = {{args.title | repr}}
figtype  = {{args.figtype | repr}}
passed   = {{args.passed | repr}}
region   = {{args.region | repr}}
regfile  = {{args.regfile | repr}}
macro    = {{args.macro | repr}}
ggs      = {{args.ggs | repr}}
devpars  = {{args.devpars | repr}}
vcfstats = {{args.vcfstats | repr}}

shell.load_config(vcfstats = vcfstats)

shell.fg.vcfstats(
	vcf      = infile,
	outdir   = outdir,
	formula  = formula,
	title    = title,
	loglevel = 'DEBUG',
	figtype  = figtype or False,
	region   = region or False,
	Region   = regfile or False,
	passed   = passed,
	macro    = macro or False,
	ggs      = ggs or False,
	devpars  = devpars or False,
	config   = config or False
)
