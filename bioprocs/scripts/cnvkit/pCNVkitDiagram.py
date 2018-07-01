from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
cnrfile  = {{in.cnrfile | quote}}
cnsfile  = {{in.cnsfile | quote}}
outfile  = {{out.outfile | quote}}
params   = {{args.params}}
nthread  = {{args.nthread | repr}}

openblas_nthr = "export OPENBLAS_NUM_THREADS={nthread}; export OMP_NUM_THREADS={nthread}; export NUMEXPR_NUM_THREADS={nthread}; export MKL_NUM_THREADS={nthread}".format(nthread = nthread)

cmd      = '{openblas}; {cnvkit} diagram \'{cnrfile}\' {params}'
params.o = outfile
if cnsfile:
	params.s = cnsfile

runcmd(cmd.format(**Box(
	openblas = openblas_nthr,
	cnvkit  = cnvkit,
	params  = cmdargs(params, equal = ' '),
	cnrfile = cnrfile
)))