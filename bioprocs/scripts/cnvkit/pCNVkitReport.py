from os import path
from pyppl import Box
from shutil import copyfile
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
cnrfile  = {{i.cnrfile | quote}}
cnsfile  = {{i.cnsfile | quote}}
stem     = {{i.cnrfile | bn | quote}}
outdir   = {{o.outdir | quote}}
nthread  = {{args.nthread | repr}}

# also report cnr, cns files
copyfile(cnrfile, path.join(outdir, path.basename(cnrfile)))
copyfile(cnsfile, path.join(outdir, path.basename(cnsfile)))

params   = {{args.params}}

openblas_nthr = "export OPENBLAS_NUM_THREADS={nthread}; export OMP_NUM_THREADS={nthread}; export NUMEXPR_NUM_THREADS={nthread}; export MKL_NUM_THREADS={nthread}".format(nthread = nthread)

if params.breaks:
	if params.breaks is True:
		params.breaks = Box()
	params.breaks.o = path.join(outdir, stem + '.breaks.txt')
	cmd = '{openblas}; {cnvkit} breaks \'{cnrfile}\' \'{cnsfile}\' {params}'
	runcmd(cmd.format(**Box(
		openblas = openblas_nthr,
		cnvkit   = cnvkit,
		params   = cmdargs(params.breaks, equal = ' '),
		cnrfile  = cnrfile,
		cnsfile  = cnsfile
	)))

if params.gainloss:
	if params.gainloss is True:
		params.gainloss = Box()
	params.gainloss.o = path.join(outdir, stem + '.gainloss.txt')
	params.gainloss.s = cnsfile
	cmd = '{openblas}; {cnvkit} gainloss \'{cnrfile}\' {params}'
	runcmd(cmd.format(**Box(
		openblas = openblas_nthr,
		cnvkit   = cnvkit,
		params   = cmdargs(params.gainloss, equal = ' '),
		cnrfile  = cnrfile,
	)))

if params.metrics:
	if params.metrics is True:
		params.metrics = Box()
	params.metrics.o = path.join(outdir, stem + '.metrics.txt')
	params.metrics.s = cnsfile
	cmd = '{openblas}; {cnvkit} metrics \'{cnrfile}\' {params}'
	runcmd(cmd.format(**Box(
		openblas = openblas_nthr,
		cnvkit   = cnvkit,
		params   = cmdargs(params.metrics, equal = ' '),
		cnrfile  = cnrfile,
	)))

if params.segmetrics:
	if params.segmetrics is True:
		params.segmetrics = Box()
	params.segmetrics.o = path.join(outdir, stem + '.segmetrics.txt')
	params.segmetrics.s = cnsfile
	cmd = '{openblas}; {cnvkit} segmetrics \'{cnrfile}\' {params}'
	runcmd(cmd.format(**Box(
		openblas = openblas_nthr,
		cnvkit   = cnvkit,
		params   = cmdargs(params.segmetrics, equal = ' '),
		cnrfile  = cnrfile,
	)))


