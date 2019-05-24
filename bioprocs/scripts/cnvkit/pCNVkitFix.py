from os import path
from pyppl import Box
from bioprocs.utils import shell2 as shell

cnvkit   = {{args.cnvkit | quote}}
tgfile   = {{i.tgfile | quote}}
atgfile  = {{i.atgfile | quote}}
rcfile   = {{i.rcfile | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params | repr}}
nthread  = {{args.nthread | repr}}
sample   = {{i.tgfile | fn | quote}}

shell.load_config(cnvkit = dict(
	_exe = cnvkit,
	_env = dict(
		OPENBLAS_NUM_THREADS = str(nthread),
		OMP_NUM_THREADS      = str(nthread),
		NUMEXPR_NUM_THREADS  = str(nthread),
		MKL_NUM_THREADS      = str(nthread)
	),
	_cwd = path.dirname(outfile)
))

params.o = outfile
params.i = sample # sample id
shell.fg.cnvkit.fix(tgfile, atgfile, rcfile, **params)
