from os import path
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

infiles = {{i.infiles | repr}}
exbaits = {{args.exbaits | repr}}
accfile = {{args.accfile | repr}}
ref     = {{args.ref | repr}}
params  = {{args.params | repr}}
prefix  = {{i.infiles | fs2name | quote}}
outdir  = {{job.outdir | quote}}
cnvkit  = {{args.cnvkit | quote}}
nthread = {{args.nthread | repr}}

# generate target file
cmd        = '{cnvkit} target {baitfile} {params}'
params_t   = params.target
params_t.o = path.join(outdir, prefix + '.bed')
runcmd(cmd.format(
	cnvkit   = cnvkit,
	baitfile = repr(exbaits),
	params   = cmdargs(params_t, equal = ' ')
))

# generate access file
if not accfile:
	accfile = path.join(outdir, prefix + '.access.bed')
	cmd = '{cnvkit} access {ref} {params}'
	params_a = params.access
	params_a.o = accfile
	runcmd(cmd.format(
		cnvkit = cnvkit,
		ref    = repr(ref),
		params = cmdargs(params_a, equal = ' ')
	))

# autobin
openblas_nthr = "export OPENBLAS_NUM_THREADS={nthread}; export OMP_NUM_THREADS={nthread}; export NUMEXPR_NUM_THREADS={nthread}; export MKL_NUM_THREADS={nthread}".format(nthread = nthread)
cmd = 'cd {outdir}; {nthr}; {cnvkit} autobin {bams} {params}'
params_b = params.autobin
params_b.t = params_t.o
params_b.g = accfile
runcmd(cmd.format(
	outdir = str(repr(outdir)),
	cnvkit = cnvkit,
	nthr   = openblas_nthr,
	bams   = ' '.join([str(repr(infile)) for infile in infiles]),
	params = cmdargs(params_b, equal = ' ')
))
