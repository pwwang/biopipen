from copy import deepcopy
from pyppl import Box
from bioprocs.utils import shell

cnvkit   = {{args.cnvkit | quote}}
infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
antifile = {{o.antifile | quote}}
target   = {{i.tgfile | quote}}
atarget  = {{i.atgfile | quote}}
nthread  = {{args.nthread | repr}}
params   = {{args.params}}

shell.TOOLS['cnvkit'] = cnvkit
envs = dict(
	OPENBLAS_NUM_THREADS = nthread,
	OMP_NUM_THREADS      = nthread,
	NUMEXPR_NUM_THREADS  = nthread,
	MKL_NUM_THREADS      = nthread
)
ckshell = shell.Shell(subcmd = True, equal = ' ', envs = envs).cnvkit

params.p  = nthread
params1   = deepcopy(params)
params2   = deepcopy(params)
params1.o = outfile
params2.o = antifile

ckshell.coverage(infile, target,  **params1).run()
ckshell.coverage(infile, atarget, **params2).run()
