from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
infile   = {{i.cnsfile | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params}}

params.o = outfile
cmd      = '{cnvkit} export vcf \'{infile}\' {params}'

runcmd(cmd.format(**Box(
	cnvkit   = cnvkit,
	params   = cmdargs(params, equal = ' '),
	infile   = infile
)))