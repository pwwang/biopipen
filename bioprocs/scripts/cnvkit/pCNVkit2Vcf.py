from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit  = {{args.cnvkit | quote}}
infile  = {{i.cnsfile | quote}}
outfile = {{o.outfile | quote}}
params  = {{args.params}}
nthread = {{args.nthread | repr}}

params.o = outfile
openblas_nthr = "export OPENBLAS_NUM_THREADS={nthread}; export OMP_NUM_THREADS={nthread}; export NUMEXPR_NUM_THREADS={nthread}; export MKL_NUM_THREADS={nthread}".format(nthread = nthread)
cmd      = '{nthr}; {cnvkit} export vcf \'{infile}\' {params}'

runcmd(cmd.format(**Box(
	nthr     = openblas_nthr,
	cnvkit   = cnvkit,
	params   = cmdargs(params, equal = ' '),
	infile   = infile
)))