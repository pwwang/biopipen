from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
nthread  = {{args.nthread | repr}}
params   = {{args.params}}

openblas_nthr = "export OPENBLAS_NUM_THREADS={nthread}; export OMP_NUM_THREADS={nthread}; export NUMEXPR_NUM_THREADS={nthread}; export MKL_NUM_THREADS={nthread}".format(nthread = nthread)

params.o = outfile
params.p = nthread
cmd      = '{openblas}; {cnvkit} segment \'{infile}\' {params}'

runcmd(cmd.format(**Box(
	openblas = openblas_nthr,
	cnvkit   = cnvkit,
	params   = cmdargs(params, equal = ' '),
	infile   = infile
)))