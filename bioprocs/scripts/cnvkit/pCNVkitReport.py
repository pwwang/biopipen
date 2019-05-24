from os import path
from pyppl import Box
from shutil import copyfile
from bioprocs.utils import shell2 as shell

cnvkit   = {{args.cnvkit | quote}}
cnrfile  = {{i.cnrfile | quote}}
cnsfile  = {{i.cnsfile | quote}}
stem     = {{i.cnrfile | bn | quote}}
outdir   = {{o.outdir | quote}}
nthread  = {{args.nthread | repr}}
params   = {{args.params | repr}}

shell.load_config(cnvkit = dict(
	_exe = cnvkit,
	_env = dict(
		OPENBLAS_NUM_THREADS = str(nthread),
		OMP_NUM_THREADS      = str(nthread),
		NUMEXPR_NUM_THREADS  = str(nthread),
		MKL_NUM_THREADS      = str(nthread)
	),
	_cwd = outdir
))

# also report cnr, cns files
shell.cp(cnrfile, path.join(outdir, path.basename(cnrfile)))
shell.cp(cnsfile, path.join(outdir, path.basename(cnsfile)))

if params.breaks:
	if params.breaks is True:
		params.breaks = Box()
	params.breaks.o = path.join(outdir, stem + '.breaks.txt')
	shell.fg.cnvkit.breaks(cnrfile, cnsfile, **params.breaks)

if params.gainloss:
	if params.gainloss is True:
		params.gainloss = Box()
	params.gainloss.o = path.join(outdir, stem + '.gainloss.txt')
	params.gainloss.s = cnsfile
	shell.fg.cnvkit.gainloss(cnrfile, **params.gainloss)

if params.metrics:
	if params.metrics is True:
		params.metrics = Box()
	params.metrics.o = path.join(outdir, stem + '.metrics.txt')
	params.metrics.s = cnsfile
	shell.fg.cnvkit.metrics(cnrfile, **params.metrics)

if params.segmetrics:
	if params.segmetrics is True:
		params.segmetrics = Box()
	params.segmetrics.o = path.join(outdir, stem + '.segmetrics.txt')
	params.segmetrics.s = cnsfile
	shell.fg.cnvkit.segmetrics(cnrfile, **params.segmetrics).run()

