from copy import deepcopy
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
infile   = {{in.infile | quote}}
outfile  = {{out.outfile | quote}}
antifile = {{out.antifile | quote}}
target   = {{in.tgfile | quote}}
atarget  = {{in.atgfile | quote}}
nthread  = {{args.nthread | repr}}
params   = {{args.params}}

openblas_nthr = "export OPENBLAS_NUM_THREADS={nthread}; export OMP_NUM_THREADS={nthread}; export NUMEXPR_NUM_THREADS={nthread}; export MKL_NUM_THREADS={nthread}".format(nthread = nthread)

params.p  = nthread
params1   = deepcopy(params)
params2   = deepcopy(params)
params1.o = outfile
params2.o = antifile
cmd       = '{openblas}; {cnvkit} coverage {infile} {target} {params1}; {cnvkit} coverage {infile} {atarget} {params2}'

runcmd(cmd.format(**Box(
	openblas = openblas_nthr,
	cnvkit   = cnvkit,
	infile   = str(repr(infile)),
	target   = str(repr(target)),
	atarget  = str(repr(atarget)),
	params1  = cmdargs(params1, equal = ' '),
	params2  = cmdargs(params2, equal = ' ')
)))