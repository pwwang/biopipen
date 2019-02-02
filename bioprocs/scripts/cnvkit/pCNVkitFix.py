from pyppl import Box
from bioprocs.utils import shell

cnvkit   = {{args.cnvkit | quote}}
tgfile   = {{i.tgfile | quote}}
atgfile  = {{i.atgfile | quote}}
rcfile   = {{i.rcfile | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params}}
nthread  = {{args.nthread | repr}}
sample   = {{i.tgfile | fn | quote}}

shell.TOOLS['cnvkit'] = cnvkit
envs = dict(
	OPENBLAS_NUM_THREADS = nthread,
	OMP_NUM_THREADS      = nthread,
	NUMEXPR_NUM_THREADS  = nthread,
	MKL_NUM_THREADS      = nthread
)
ckshell = shell.Shell(subcmd = True, equal = ' ', envs = envs).cnvkit


params.o = outfile
params.i = sample # sample id
ckshell.fix(tgfile, atgfile, rcfile, **params).run()
