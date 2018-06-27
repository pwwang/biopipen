from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
infile   = {{in.infile | quote}}
rcfile   = {{in.rcfile | quote}}
outfile  = {{out.outfile | quote}}
params   = {{args.params}}

mtfile   = '{{job.outdir}}/cnvkit_mt'
open(mtfile, 'w').close()

params.o = outfile
cmd      = '{cnvkit} fix \'{infile}\' \'{mtfile}\' \'{rcfile}\' {params}'

runcmd(cmd.format(**Box(
	cnvkit   = cnvkit,
	mtfile   = mtfile,
	params   = cmdargs(params, equal = ' '),
	infile   = infile,
	rcfile   = rcfile
)))