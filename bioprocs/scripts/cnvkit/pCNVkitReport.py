from os import path
from pyppl import Box
from shutil import copyfile
from bioprocs.utils import shell

cnvkit   = {{args.cnvkit | quote}}
cnrfile  = {{i.cnrfile | quote}}
cnsfile  = {{i.cnsfile | quote}}
stem     = {{i.cnrfile | bn | quote}}
outdir   = {{o.outdir | quote}}
nthread  = {{args.nthread | repr}}
params   = {{args.params}}

shell.TOOLS['cnvkit'] = cnvkit
envs = dict(
	OPENBLAS_NUM_THREADS = nthread,
	OMP_NUM_THREADS      = nthread,
	NUMEXPR_NUM_THREADS  = nthread,
	MKL_NUM_THREADS      = nthread
)
ckshell = shell.Shell(subcmd = True, equal = ' ', envs = envs, cwd = outdir).cnvkit

# also report cnr, cns files
shell.cp(cnrfile, path.join(outdir, path.basename(cnrfile)))
shell.cp(cnsfile, path.join(outdir, path.basename(cnsfile)))

if params.breaks:
	if params.breaks is True:
		params.breaks = Box()
	params.breaks.o = path.join(outdir, stem + '.breaks.txt')
	ckshell.breaks(cnrfile, cnsfile, **params.breaks).run()

if params.gainloss:
	if params.gainloss is True:
		params.gainloss = Box()
	params.gainloss.o = path.join(outdir, stem + '.gainloss.txt')
	params.gainloss.s = cnsfile
	ckshell.gainloss(cnrfile, **params.gainloss).run()

if params.metrics:
	if params.metrics is True:
		params.metrics = Box()
	params.metrics.o = path.join(outdir, stem + '.metrics.txt')
	params.metrics.s = cnsfile
	ckshell.metrics(cnrfile, **params.metrics).run()

if params.segmetrics:
	if params.segmetrics is True:
		params.segmetrics = Box()
	params.segmetrics.o = path.join(outdir, stem + '.segmetrics.txt')
	params.segmetrics.s = cnsfile
	ckshell.segmetrics(cnrfile, **params.segmetrics).run()

