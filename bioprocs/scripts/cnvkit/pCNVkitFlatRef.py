from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
infile   = {{in.tgfile | quote}}
atgfile  = {{in.atgfile | quote}}
ref      = {{args.ref | quote}}
outfile  = {{out.outfile | quote}}
params   = {{args.params}}

params.o = outfile
params.f = ref
params.t = infile
params.a = atgfile
cmd      = '{cnvkit} reference {params}'

runcmd(cmd.format(**Box(
	cnvkit   = cnvkit,
	params   = cmdargs(params, equal = ' ')
)))