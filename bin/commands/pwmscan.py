#!/usr/bin/env python
from pyppl import PyPPL
from bioprocs import params
from bioprocs.gene import pGeneNameNorm
from bioaggrs.tfbs import aTfbsTfPC, aTfbsTfRC, aTfbsTfP, aTfbsTfR, aTfbsPC, aTfbsRC, aTfbsP, aTfbsR

#params.prefix('-')

params.tflist.show     = True
params.tfmotifs.show   = True
params.consvdir.show   = True

params.input           = "tf"
params.input.desc      = "The input content. tf or motif."

params.infile.required = True
params.infile.desc     = "The input file. Whether TF or motif list."

params.target          = 'gene'
params.target.desc     = "What's the target. gene or region"

params.tfile.required  = True
params.tfile.desc      = 'The target file. If target is region, using BED format.'

params.pval            = 1e-4
params.pval.desc       = "The p-value cutoff for PWM scan."

params.consvp          = 5e-2
params.consvp.desc     = "The conservation p-value cutoff. 0 to disable conservation report."

params.outdir.required = True
params.outdir.desc     = "The output directory."

params.tfnorm          = True
params.tfnorm.desc     = 'Whether normalize TF names.'

params.nthread         = 1
params.nthread.desc    = 'Number of threads used for scanning.'

params = params.parse().toDict()

# check params
if not params.input in ['tf', 'motif']:
	raise ValueError('Unexpected "input" type, should be one of ["tf", "motif"]')

if not params.target in ['gene', 'region']:
	raise ValueError('Unexpected "target" type, should be one of ["gene", "region"]')

if params.tfnorm and params.input == 'tf':
	pGeneNameNorm.input = [params.infile]
	PyPPL().start(pGeneNameNorm).run()
	params.infile = pGeneNameNorm.channel.get()

if params.input == 'tf' and params.target == 'gene' and params.consvp > 0:
	aggr = aTfbsTfPC
elif params.input == 'tf' and params.target == 'gene' and params.consvp == 0:
	aggr = aTfbsTfP
elif params.input == 'tf' and params.target == 'region' and params.consvp > 0:
	aggr = aTfbsTfRC
elif params.input == 'tf' and params.target == 'region' and params.consvp == 0:
	aggr = aTfbsTfR
if params.input == 'motif' and params.target == 'gene' and params.consvp > 0:
	aggr = aTfbsPC
elif params.input == 'motif' and params.target == 'gene' and params.consvp == 0:
	aggr = aTfbsP
elif params.input == 'motif' and params.target == 'region' and params.consvp > 0:
	aggr = aTfbsRC
elif params.input == 'motif' and params.target == 'region' and params.consvp == 0:
	aggr = aTfbsR

if params.input == 'tf':
	aggr.input = [[params.infile], [params.tfile], [params.tflist]]
else:
	aggr.input = [[params.infile], [params.tfile]]

if params.consvp > 0:
	aggr.args.cpval = params.consvp
	aggr.pConsvPerm.args.consvdir = params.consvdir

aggr.args.pval               = params.pval
aggr.args.tfmotifs           = params.tfmotifs
aggr.cclean                  = True
aggr.exdir                   = params.outdir
aggr.pMotifScan.args.nthread = params.nthread

PyPPL().start(aggr).run()







