from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
cnrfile  = {{in.cnrfile | quote}}
cnsfile  = {{in.cnsfile | quote}}
outfile  = {{out.outfile | quote}}
params   = {{args.params}}
cmd      = '{cnvkit} diagram \'{cnrfile}\' {params}'

params.o = outfile
if cnsfile:
	params.s = cnsfile

runcmd(cmd.format(**Box(
	cnvkit  = cnvkit,
	params  = cmdargs(params, equal = ' '),
	cnrfile = cnrfile
)))