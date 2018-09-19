from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
infile   = {{i.tgfile | quote}}
atgfile  = {{i.atgfile | quote}}
ref      = {{args.ref | quote}}
outfile  = {{o.outfile | quote}}
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