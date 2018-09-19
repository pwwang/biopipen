from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
infiles  = {{i.infiles | asquote | quote}}
ref      = {{args.ref | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params}}

params.o = outfile
params.f = ref
cmd      = '{cnvkit} reference {infiles} {params}'

runcmd(cmd.format(**Box(
	cnvkit   = cnvkit,
	params   = cmdargs(params, equal = ' '),
	infiles  = infiles
)))