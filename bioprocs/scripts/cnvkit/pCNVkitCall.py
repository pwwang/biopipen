from pyppl import Box
from bioprocs.utils import shell

cnvkit   = {{args.cnvkit | quote}}
infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params}}

shell.TOOLS['cnvkit'] = cnvkit
# envs = dict(
# 	OPENBLAS_NUM_THREADS = nthread,
# 	OMP_NUM_THREADS      = nthread,
# 	NUMEXPR_NUM_THREADS  = nthread,
# 	MKL_NUM_THREADS      = nthread
# )
ckshell = shell.Shell(subcmd = True, equal = ' ').cnvkit

params.o = outfile
ckshell.call(infile, **params).run()
