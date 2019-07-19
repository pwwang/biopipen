#!/usr/bin/env python
"""Do PWM scan"""
from pyppl import PyPPL, Box
from bioprocs import params

params._prefix = '-'

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

def main():
	opts = params._parse(dict_wrapper = Box)

	from bioprocs.gene import pGeneNameNorm
	from bioprocs.sets.tfbs import aTfbsTfPC, aTfbsTfRC, aTfbsTfP, aTfbsTfR, \
		aTfbsPC, aTfbsRC, aTfbsP, aTfbsR

	# check params
	if not opts.input in ['tf', 'motif']:
		raise ValueError('Unexpected "input" type, should be one of ["tf", "motif"]')

	if not opts.target in ['gene', 'region']:
		raise ValueError('Unexpected "target" type, should be one of ["gene", "region"]')

	if opts.tfnorm and opts.input == 'tf':
		pGeneNameNorm.input = [opts.infile]
		PyPPL().start(pGeneNameNorm).run()
		opts.infile = pGeneNameNorm.channel.get()

	if opts.input == 'tf' and opts.target == 'gene' and opts.consvp > 0:
		aggr = aTfbsTfPC
	elif opts.input == 'tf' and opts.target == 'gene' and opts.consvp == 0:
		aggr = aTfbsTfP
	elif opts.input == 'tf' and opts.target == 'region' and opts.consvp > 0:
		aggr = aTfbsTfRC
	elif opts.input == 'tf' and opts.target == 'region' and opts.consvp == 0:
		aggr = aTfbsTfR
	if opts.input == 'motif' and opts.target == 'gene' and opts.consvp > 0:
		aggr = aTfbsPC
	elif opts.input == 'motif' and opts.target == 'gene' and opts.consvp == 0:
		aggr = aTfbsP
	elif opts.input == 'motif' and opts.target == 'region' and opts.consvp > 0:
		aggr = aTfbsRC
	elif opts.input == 'motif' and opts.target == 'region' and opts.consvp == 0:
		aggr = aTfbsR

	if opts.input == 'tf':
		aggr.input = [[opts.infile], [opts.tfile], [opts.tflist]]
	else:
		aggr.input = [[opts.infile], [opts.tfile]]

	if opts.consvp > 0:
		aggr.args.cpval = opts.consvp
		aggr.pConsvPerm.args.consvdir = opts.consvdir

	aggr.args.pval               = opts.pval
	aggr.args.tfmotifs           = opts.tfmotifs
	aggr.acache                  = True
	aggr.exdir                   = opts.outdir
	aggr.pMotifScan.args.nthread = opts.nthread

	PyPPL().start(aggr).run()

if __name__ == "__main__":
	main()





