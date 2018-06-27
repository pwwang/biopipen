from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
#index    = {{in.index | quote}}
ref      = {{args.ref | quote}}
outfile  = {{out.outfile | quote}}
params   = {{args.params}}
params.o = outfile
cmd      = '{cnvkit} access \'{infile}\' {params}'

runcmd(cmd.format(**Box(
	cnvkit  = cnvkit,
	params  = cmdargs(params, equal = ' '),
	infile  = ref
)))