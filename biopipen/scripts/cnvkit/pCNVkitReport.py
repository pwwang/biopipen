from os import path
from diot import Diot
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
		params.breaks = Diot()
	params.breaks.o = path.join(outdir, stem + '.breaks.txt')
	shell.cnvkit.breaks(cnrfile, cnsfile, **params.breaks).fg

if params.gainloss:
	if params.gainloss is True:
		params.gainloss = Diot()
	params.gainloss.o = path.join(outdir, stem + '.gainloss.txt')
	params.gainloss.s = cnsfile
	shell.cnvkit.gainloss(cnrfile, **params.gainloss).fg

if params.metrics:
	if params.metrics is True:
		params.metrics = Diot()
	params.metrics.o = path.join(outdir, stem + '.metrics.txt')
	params.metrics.s = cnsfile
	shell.cnvkit.metrics(cnrfile, **params.metrics).fg

if params.segmetrics:
	if params.segmetrics is True:
		params.segmetrics = Diot()
	params.segmetrics.o = path.join(outdir, stem + '.segmetrics.txt')
	params.segmetrics.s = cnsfile
	shell.cnvkit.segmetrics(cnrfile, **params.segmetrics).fg
