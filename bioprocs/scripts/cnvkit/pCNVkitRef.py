from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
infiles  = {{in.infiles | asquote | quote}}
ref      = {{args.ref | quote}}
outfile  = {{out.outfile | quote}}
params   = {{args.params}}

params.o = outfile
params.f = ref
cmd      = '{cnvkit} reference {infiles} {params}'

runcmd(cmd.format(**Box(
	cnvkit   = cnvkit,
	params   = cmdargs(params, equal = ' '),
	infiles  = infiles
)))