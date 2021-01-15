#!/usr/bin/env python
"""Call mutatoins from 2nd-gen sequencing data"""
# pylint: disable=invalid-name
from os import path
from pathlib import Path
from pyppl import PyPPL, Channel
from pyppl.logger import logger
from pyparam import Params
from biopipen import params as bp_params

HERE = Path(__file__).parent.resolve()
ARGS_FILE = HERE / 'mutcall.args.toml'

help_group = 'SCRIPTS'
params = Params(desc=__doc__)

params.from_file(ARGS_FILE)


def main(opts): # pylint: disable=too-many-statements,assigning-non-slot
    """The main entry point of the pipeline"""
    opts = bp_params.parse() | opts

    from biopipen.common import pFiles2Dir
    from biopipen.sambam import pBam2Gmut, pBamPair2Smut
    from biopipen.utils.sampleinfo import SampleInfo2 as SampleInfo
    from biopipen.sets.wxs import aPrepareBam, aBam2SCNV, aBam2GCNV

    starts = []
    saminfo = SampleInfo(opts.saminfo)
    aPrepareBam.pFastq2Sam.args.tool = opts.aligner
    aPrepareBam.pFastqTrim.args.tool = opts.trimmer
    aPrepareBam.args.nthread = opts.nthread

    if (aPrepareBam.pSam2Bam.args.tool == 'elprep' and
            aPrepareBam.pSam2Bam.args.steps.recal):
        aPrepareBam.modules.norecal()
    if opts.compress:
        aPrepareBam.args.gz = True
        aPrepareBam.pFastq2Sam.args.outfmt = 'bam'

    pBamDir = pFiles2Dir # pylint: disable=invalid-name
    pBamDir.runner = 'local'
    if opts.intype == 'ebam':
        logger.title(
            'Mutcall: External BAMs received, try to refactorize them.'
        )
        #aPrepareBam.input = [
        #   Channel.fromPattern(path.join(opts.indir, '*.bam'))]
        aPrepareBam.modules.ebam(restore=False)
        aPrepareBam.input = Channel.create(saminfo.toChannel(
            opts.indir)).unique()

        pBamDir.depends = aPrepareBam
        pBamDir.input = lambda ch: [ch.flatten()]

        starts.append(aPrepareBam)

    elif opts.intype in ('fq', 'fastq'):
        # pair-end fastq files
        # *.fq, *.fq.gz *.fastq, *.fastq.gz
        # sample info should be:
        # +--------------+----------+---------+
        # | Sample       | Patient  | Group   |
        # | x_Tumor.bam  | x        | TUMOR   |
        # | x_Normal.bam | x        | NORMAL  |
        # | ...          | ...      | ...     |
        # +--------------+----------+---------+
        # corresponding fastq files would be:
        # x_Tumor_1.fq(.gz)  / x_Tumor_1.fastq(.gz)
        # x_Tumor_2.fq(.gz)  / x_Tumor_2.fastq(.gz)
        # x_Normal_1.fq(.gz) / x_Normal_1.fastq(.gz)
        # x_Normal_2.fq(.gz) / x_Normal_2.fastq(.gz)
        def bam2fqpair(fastq):
            fqdir = path.dirname(fastq)
            bname = path.splitext(path.basename(fastq))[0]
            exts1 = ['_1.fq', '_1.fq.gz', '_1.fq.gz', '_1.fastq.gz']
            exts2 = ['_2.fq', '_2.fq.gz', '_2.fq.gz', '_2.fastq.gz']
            fqfiles1 = [path.join(fqdir, bname + ext) for ext in exts1]
            fqfile1 = [fqfile for fqfile in fqfiles1 if path.isfile(fqfile)][0]
            fqfiles2 = [path.join(fqdir, bname + ext) for ext in exts2]
            fqfile2 = [fqfile for fqfile in fqfiles2 if path.isfile(fqfile)][0]
            return fqfile1, fqfile2

        aPrepareBam.modules.fastq()
        aPrepareBam.input = [
            bam2fqpair(fastq) for fastq in saminfo.toChannel(opts.indir)
        ]
        pBamDir.depends = (aPrepareBam.pSam2Bam
                           if aPrepareBam.pSam2Bam.args.tool == 'elprep'
                           else aPrepareBam.pBamRecal)
        pBamDir.input = lambda ch: [ch.flatten()]

        starts.append(aPrepareBam)
    else:
        pBamDir.input = [saminfo.toChannel(opts.indir)]
        starts.append(pBamDir)

    if 'germ' in opts.muts:
        pBam2Gmut.depends = pBamDir
        pBam2Gmut.input = lambda ch: ch.expand(0, "*.bam")
        pBam2Gmut.args.nthread = opts.nthread
        pBam2Gmut.config.export_dir = path.join(opts.outdir, 'germline')
    if 'soma' in opts.muts:
        pBamPair2Smut.depends = pBamDir
        pBamPair2Smut.input = lambda ch: saminfo.toChannel(ch.get(),
                                                           paired=True)
        pBamPair2Smut.args.nthread = opts.nthread
        pBamPair2Smut.config.export_dir = path.join(opts.outdir, 'somatic')
    if 'scnv' in opts.muts:
        aBam2SCNV.modules.plots()
        aBam2SCNV.pBamDir.depends = pBamDir
        aBam2SCNV.pSampleInfo.input = [opts.saminfo]
        aBam2SCNV.args.nthread = opts.nthread
        aBam2SCNV.ends.config.export_dir = path.join(opts.outdir, 'scnv')
        starts.append(aBam2SCNV)
    if 'gcnv' in opts.muts:
        aBam2GCNV.modules.plots()
        aBam2GCNV.pBamDir.depends = pBamDir
        aBam2GCNV.pSampleInfo.input = [opts.saminfo]
        aBam2GCNV.args.nthread = opts.nthread
        aBam2GCNV.ends.config.export_dir = path.join(opts.outdir, 'gcnv')
        starts.append(aBam2GCNV)

    ppl = PyPPL(forks=int(opts.forks),
                ppldir=opts.ppldir,
                cache='force' if opts.forcecache else True,
                errhow=opts.errhow,
                logger_level=opts.loglevel,
                logger_file=opts.logfile).start(starts)
    if opts.flowchart:
        ppl.flowchart(fcfile=opts.flowchart)
    ppl.run(opts.runner)

if __name__ == "__main__":
    main(params.parse())
