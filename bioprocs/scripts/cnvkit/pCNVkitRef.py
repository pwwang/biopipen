from pyppl import Box
from bioprocs.utils import shell

cnvkit   = {{args.cnvkit | quote}}
infiles  = {{i.infiles | repr}}
ref      = {{args.ref | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params}}
nthread  = {{args.nthread | repr}}

shell.TOOLS['cnvkit'] = cnvkit
envs = dict(
	OPENBLAS_NUM_THREADS = nthread,
	OMP_NUM_THREADS      = nthread,
	NUMEXPR_NUM_THREADS  = nthread,
	MKL_NUM_THREADS      = nthread
)
ckshell = shell.Shell(subcmd = True, equal = ' ', envs = envs).cnvkit

params.o = outfile
params.f = ref
ckshell.reference(*infiles, **params).run()
