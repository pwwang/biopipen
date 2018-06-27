from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
infile   = {{in.acfile | quote}}
outfile  = {{out.outfile | quote}}
params   = {{args.params}}
params.o = outfile
cmd      = '{cnvkit} target \'{infile}\' {params}'

runcmd(cmd.format(**Box(
	cnvkit  = cnvkit,
	params  = cmdargs(params, equal = ' '),
	infile  = infile
)))