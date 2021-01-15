#!/usr/bin/env python
"""Do PWM scan"""
from pyppl import PyPPL
from pyparam import Params
from biopipen import params as bp_params

# pylint: disable=invalid-name
help_group = 'SCRIPTS'
params = Params(desc=__doc__)

params.add_param(bp_params.get_param('tflist'), show=True)
params.add_param(bp_params.get_param('tfmotifs'), show=True)
params.add_param(bp_params.get_param('consvdir'), show=True)

params.add_param("input", default="tf", required=True, force=True,
                 desc="The input content. tf or motif.")
params.add_param("target", default='gene',
                 desc="What's the target. gene or region")
params.add_param("tfile", required=True, type='path',
                 desc='The target file. If target is region, using BED format.')
params.add_param("pval", default=1e-4, desc="The p-value cutoff for PWM scan.")
params.add_param("consvp", default=5e-2,
                 desc=("The conservation p-value cutoff. "
                       "0 to disable conservation report."))
params.add_param("outdir", required=True, type="path", force=True,
                 desc="The output directory.")
params.add_param("tfnorm", default=True, desc='Whether normalize TF names.')
params.add_param("nthread", default=1, force=True,
                 desc='Number of threads used for scanning.')

def main(opts): # pylint: disable=too-many-branches
    """Main function"""
    opts = bp_params.parse() | opts

    from biopipen.gene import pGeneNameNorm
    from biopipen.sets.tfbs import (aTfbsTfPC,
                                    aTfbsTfRC,
                                    aTfbsTfP,
                                    aTfbsTfR,
                                    aTfbsPC,
                                    aTfbsRC,
                                    aTfbsP,
                                    aTfbsR)

    # check params
    if not opts.input in ['tf', 'motif']:
        raise ValueError(
            'Unexpected "input" type, should be one of ["tf", "motif"]')

    if not opts.target in ['gene', 'region']:
        raise ValueError(
            'Unexpected "target" type, should be one of ["gene", "region"]')

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

    aggr.args.pval = opts.pval
    aggr.args.tfmotifs = opts.tfmotifs
    aggr.acache = True
    aggr.exdir = opts.outdir
    aggr.pMotifScan.args.nthread = opts.nthread

    PyPPL().start(aggr).run()


if __name__ == "__main__":
    main(params.parse())
