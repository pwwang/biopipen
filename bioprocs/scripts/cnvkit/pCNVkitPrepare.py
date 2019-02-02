from os import path
from pyppl import Box
from bioprocs.utils import shell, runcmd
from bioprocs.utils.reference import bamIndex

infiles = {{i.infiles | repr}}
exbaits = {{args.exbaits | repr}}
accfile = {{args.accfile | repr}}
ref     = {{args.ref | repr}}
params  = {{args.params | repr}}
prefix  = {{i.infiles | fs2name | quote}}
outdir  = {{job.outdir | quote}}
cnvkit  = {{args.cnvkit | quote}}
nthread = {{args.nthread | repr}}

for infile in infiles:
	bamIndex(infile)

shell.TOOLS['cnvkit'] = cnvkit
envs = dict(
	OPENBLAS_NUM_THREADS = nthread,
	OMP_NUM_THREADS      = nthread,
	NUMEXPR_NUM_THREADS  = nthread,
	MKL_NUM_THREADS      = nthread
)
ckshell = shell.Shell(subcmd = True, equal = ' ', envs = envs, cwd = outdir).cnvkit

# generate target file
params_t   = params.target
params_t.o = path.join(outdir, prefix + '.bed')
ckshell.target(exbaits, **params_t).run()

# generate access file
if not accfile:
	accfile = path.join(outdir, prefix + '.access.bed')
	params_a = params.access
	params_a.o = accfile
	ckshell.access(ref, **params_a).run()

# autobin
params_b = params.autobin
params_b.t = params_t.o
params_b.g = accfile
runcmd('cd {}; {}'.format(
	shell.shquote(outdir),
	ckshell.autobin(*infiles, **params_b).cmd
))
