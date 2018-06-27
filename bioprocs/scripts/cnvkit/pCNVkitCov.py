from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
infile   = {{in.infile | quote}}
outfile  = {{out.outfile | quote}}
target   = {{in.tgfile | quote}}
nthread  = {{args.nthread | repr}}
params   = {{args.params}}

openblas_nthr = "export OPENBLAS_NUM_THREADS={nthread}; export OMP_NUM_THREADS={nthread}; export NUMEXPR_NUM_THREADS={nthread}; export MKL_NUM_THREADS={nthread}".format(nthread = nthread)

params.o = outfile
params.p = nthread
cmd      = '{openblas}; {cnvkit} coverage \'{infile}\' \'{target}\' {params}'

runcmd(cmd.format(**Box(
	openblas = openblas_nthr,
	cnvkit   = cnvkit,
	target   = target,
	params   = cmdargs(params, equal = ' '),
	infile   = infile
)))