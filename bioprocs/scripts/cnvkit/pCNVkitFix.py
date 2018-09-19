from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
tgfile   = {{i.tgfile | quote}}
atgfile  = {{i.atgfile | quote}}
rcfile   = {{i.rcfile | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params}}
nthread  = {{args.nthread | repr}}
sample   = {{i.tgfile | fn | quote}}

openblas_nthr = "export OPENBLAS_NUM_THREADS={nthread}; export OMP_NUM_THREADS={nthread}; export NUMEXPR_NUM_THREADS={nthread}; export MKL_NUM_THREADS={nthread}".format(nthread = nthread)

params.o = outfile
params.i = sample # sample id
cmd      = '{openblas}; {cnvkit} fix {tgfile} {atgfile} {rcfile} {params}'

runcmd(cmd.format(**Box(
	openblas = openblas_nthr,
	cnvkit   = cnvkit,
	tgfile   = str(repr(tgfile)),
	atgfile  = str(repr(atgfile)),
	rcfile   = str(repr(rcfile)),
	params   = cmdargs(params, equal = ' ')
)))