from os import path
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
cnrfile  = {{in.cnrfile | quote}}
cnsfile  = {{in.cnsfile | quote}}
stem     = {{in.cnrfile | bn | quote}}
outdir   = {{out.outdir | quote}}

params   = {{args.params}}

if params.breaks:
	if params.breaks is True:
		params.breaks = Box()
	params.breaks.o = path.join(outdir, stem + '.breaks.txt')
	cmd = '{cnvkit} breaks \'{cnrfile}\' \'{cnsfile}\' {params}'
	runcmd(cmd.format(**Box(
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
	cmd = '{cnvkit} gainloss \'{cnrfile}\' {params}'
	runcmd(cmd.format(**Box(
		cnvkit   = cnvkit,
		params   = cmdargs(params.gainloss, equal = ' '),
		cnrfile  = cnrfile,
	)))

if params.metrics:
	if params.metrics is True:
		params.metrics = Box()
	params.metrics.o = path.join(outdir, stem + '.metrics.txt')
	params.metrics.s = cnsfile
	cmd = '{cnvkit} metrics \'{cnrfile}\' {params}'
	runcmd(cmd.format(**Box(
		cnvkit   = cnvkit,
		params   = cmdargs(params.metrics, equal = ' '),
		cnrfile  = cnrfile,
	)))

if params.segmetrics:
	if params.segmetrics is True:
		params.segmetrics = Box()
	params.segmetrics.o = path.join(outdir, stem + '.segmetrics.txt')
	params.segmetrics.s = cnsfile
	cmd = '{cnvkit} segmetrics \'{cnrfile}\' {params}'
	runcmd(cmd.format(**Box(
		cnvkit   = cnvkit,
		params   = cmdargs(params.segmetrics, equal = ' '),
		cnrfile  = cnrfile,
	)))


