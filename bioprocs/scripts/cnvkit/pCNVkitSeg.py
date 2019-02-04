from pyppl import Box
from bioprocs.utils import shell

cnvkit   = {{args.cnvkit | quote}}
infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
nthread  = {{args.nthread | repr}}
params   = {{args.params}}

shell.TOOLS['cnvkit'] = cnvkit
envs = dict(
	OPENBLAS_NUM_THREADS = 1,
	OMP_NUM_THREADS      = 1,
	NUMEXPR_NUM_THREADS  = 1,
	MKL_NUM_THREADS      = 1
)
ckshell = shell.Shell(subcmd = True, equal = ' ', envs = envs).cnvkit

params.o = outfile
params.p = nthread
ckshell.segment(infile, **params).run()
