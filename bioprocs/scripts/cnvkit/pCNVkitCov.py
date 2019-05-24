from os import path
from pyppl import Box
from bioprocs.utils import shell2 as shell

cnvkit   = {{args.cnvkit | quote}}
infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
antifile = {{o.antifile | quote}}
target   = {{i.tgfile | quote}}
atarget  = {{i.atgfile | quote}}
nthread  = {{args.nthread | repr}}
params   = {{args.params | repr}}

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

params.p  = nthread
params1   = params.copy()
params2   = params.copy()
params1.o = outfile
params2.o = antifile

shell.fg.cnvkit.coverage(infile, target,  **params1)
shell.fg.cnvkit.coverage(infile, atarget, **params2)
