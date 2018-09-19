from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params}}

params.o = outfile
cmd      = '{cnvkit} call \'{infile}\' {params}'

runcmd(cmd.format(**Box(
	cnvkit   = cnvkit,
	params   = cmdargs(params, equal = ' '),
	infile   = infile
)))