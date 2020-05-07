#!/usr/bin/env python
"""Run DEG analysis for expression data:
1. Statistics of given expression matrix
2. Batch effect correction if Batch is found in sample info file
3. Statistics of corrected expression matrix
4. DEG calling
5. GSEA analysis
6. Enrichr analysis
"""
# pylint: disable=invalid-name,assigning-non-slot

# Core information needed:
# 1. expression matrix, genes as rows and samples as columns,
#    must be count data.
# 2. sample information file, specifying Samples, Patients (if paired analysis)
#    and Group (comparison)

from os import path
from uuid import uuid5, NAMESPACE_URL

from pyppl import PyPPL
from diot import Diot
from bioprocs import params

params._desc = __doc__
params.exprmat.required = True
params.exprmat.desc = 'The expression matrix.'
params.meta.type = str
params.meta.desc = [
    'The meta information for genes. Could be one of:',
    '   - A matrix file with exactly the same rownames as `exprmat`,',
    '   - A column name list indicating position of the metadata in `exprmat`,',
    '   - A column indexes (1-based, not including rownames) '
    'indicating position of the metadata in `exprmat`,',
    '       - You can short-formatted indexes, i.e. "0-3,5" for 0,1,2,3,5',
    '   - If not provided, the rownames of `exprmat` must be gene symbols '
    'for enrichment analysis.'
]
params.gscol.desc = [
    'The gene symbol column, '
    'could be column name or index (1-based, not including rownames) '
    'relative to `meta`',
    '   - If not provided, will assume rownames of `exprmat` are gene symbols.'
]
params.saminfo.required = True
params.saminfo.desc = 'The sample information file.'
params.ppldir = './workdir'
params.ppldir.desc = 'The pipeline directory.'
params.exdir.desc = 'Where to export the result files.'
params.caller = 'deseq2'
params.caller.desc = 'The DEG caller'
params.runner = 'local'
params.runner.desc = 'The runner to run the processes.'
params.nthread = 1
params.nthread.desc = '# threads to use.'
params.enlibs = ['KEGG_2019_Human']
params.enlibs.desc = 'The enrichment analysis libraries'
params.enincl.desc = ('A lambda function of a TsvRecord to '
                      'include genes in Enrichr analysis')
params.cutoff = 0.05
params.cutoff.desc = 'The cutoff for DEGs'
params.cutoff.callback = lambda opt: opt.set_value({
    "by": "q",
    "value": opt.value
} if not isinstance(opt.value, dict) else opt.value)
params.filter = ""
params.filter.desc = "A filter for the expression matrix"
params.filter.callback = lambda opt: opt.set_value(opt.value or None)
params.skips = []
params.skips.desc = [
    'Skip some steps, including: stats, batch, call, gsea and enrich',
    '   - stats: Statistics for the expression profile. '
    'If batch effect correction is enabled, '
    'will do for pre- and post-correction stats.',
    '   - batch: Batch effect correction.',
    '   - call: DEG calling',
    '   - gsea: The GSEA analysis (coming) (requires call)',
    '   - enrich: Enrichr enrichment analysis (requires call)'
]
params.skips.callback = lambda opt: (
    'Unknown steps to skip.'
    if set(opt.value) - {'stats', 'batch', 'call', 'gsea', 'enrich'}
    else True
)
params.report.desc = 'The report file or a directory to save the file.'
params.title = 'Reports for DEG analysis'
params.title.desc = 'The title of the report'

def expand_numbers(nums):
    """Expand short formatted numbers like '0-3,5' to [0,1,2,3,5]"""
    parts = [part.strip() for part in nums.split(',')]
    ret = []
    for part in parts:
        if '-' not in part:
            ret.append(int(part))
        else:
            pt1, pt2 = part.split('-', 1)
            pt1, pt2 = int(pt1), int(pt2)
            ret.extend(range(pt1, pt2+1))
    return ret


def main(): # pylint: disable=too-many-locals,too-many-statements,too-many-branches
    """Main entry point"""
    opts = params._parse(dict_wrapper=Diot)
    dirindex = 1

    from bioprocs.tsv import pTsvColSelect
    from bioprocs.gsea import pEnrichr
    from bioprocs.rnaseq import pExprStats, pRNASeqDEG, pBatchEffect
    from bioprocs.common import pFile2Proc
    from bioprocs.utils.sampleinfo import SampleInfo2
    from bioprocs.utils.tsvio2 import TsvReader

    saminfo = SampleInfo2(opts.saminfo)
    samples = saminfo.all_samples(unique=True)
    reader = TsvReader(opts.exprmat)
    reader.close()

    pTsvColSelectSamples = pTsvColSelect.copy(desc='Select expression values '
                                              'for samples from the sample '
                                              'information file.')
    pTsvColSelectSamples.input = [opts.exprmat]
    pTsvColSelectSamples.args.cols = [reader.cnames[0]] + samples
    # in case we have multiple comparisons against the same expression matrix
    pTsvColSelectSamples.tag = str(uuid5(
        NAMESPACE_URL,
        ':'.join(pTsvColSelectSamples.args.cols)
    ))[:8]
    starts = [pTsvColSelectSamples]

    if opts.meta and path.isfile(opts.meta):
        pMeta = pFile2Proc.copy(desc='Get meta data file')
        pMeta.input = [opts.meta]
    else:
        pMeta = pTsvColSelect.copy(desc='Extract meta information '
                                   'from the expression matrix')
        pMeta.input = [opts.exprmat]
        pMeta.args.cols = [0] + (expand_numbers(opts.meta)
                                 if opts.meta and opts.meta[0].isdigit()
                                 else (opts.meta or []))
    starts.append(pMeta)

    if opts.filter and 'stats' in opts.skips:
        raise ValueError('Step stats cannot be skipped '
                         'for expression filtering')
    if 'stats' not in opts.skips:
        pExprStats.depends = pTsvColSelectSamples
        pExprStats.input = lambda ch: ch.cbind(opts.saminfo)
        pExprStats.args.tsform = 'function(x) log2(x+1)'
        pExprStats.args.params.histogram.bins = 100
        pExprStats.args.filter = opts.filter
        if opts.exdir:
            pExprStats.config.export_dir = path.join(
                opts.exdir,
                '%d. Expression profile statistics' % dirindex
            )
            dirindex += 1
        if opts.filter:
            pTsvColSelectSamples = pExprStats

    if saminfo.all_batches() and 'batch' not in opts.skips:
        pBatchEffect.depends = pTsvColSelectSamples
        if opts.filter:
            pBatchEffect.input = lambda ch: (ch.expand(pattern='*')
                                             .filter_col(lambda x:
                                                         not x.endswith('.png'))
                                             .cbind(opts.saminfo))
        else:
            pBatchEffect.input = lambda ch: ch.cbind(opts.saminfo)

        if 'stats' not in opts.skips:
            pExprStatsPostBatch = pExprStats.copy(
                desc='Expression profile '
                'statistics after batch effect correction'
            )
            pExprStatsPostBatch.depends = pBatchEffect
            pExprStatsPostBatch.input = lambda ch: ch.outfile.cbind(
                opts.saminfo
            )
            if opts.exdir:
                pExprStatsPostBatch.config.export_dir = path.join(
                    opts.exdir,
                    '%d. Expression profile statistics '
                    'after batch effect correction' % dirindex
                )
                dirindex += 1

        pTsvColSelectSamples = pBatchEffect

    if 'call' not in opts.skips:
        pRNASeqDEG.depends = pTsvColSelectSamples, pMeta
        pRNASeqDEG.input = (
            lambda ch1, ch2: ch1.outfile.cbind(opts.saminfo, ch2)
            if hasattr(ch1, 'outfile')
            # make sure it is .txt
            else ch1.expand(
                pattern="*%s" % path.splitext(opts.exprmat)[1]
            ).cbind(opts.saminfo, ch2))
        pRNASeqDEG.args.tool = opts.caller
        pRNASeqDEG.args.gscol = opts.gscol
        pRNASeqDEG.args.meta = opts.meta
        pRNASeqDEG.args.cutoff = opts.cutoff
        if opts.exdir:
            pRNASeqDEG.config.export_dir = path.join(
                opts.exdir,
                '%d. Differentially-expression genes' % dirindex
            )
            dirindex += 1

    if 'gsea' not in opts.skips:
        # to be implemented
        pass

    if 'enrich' not in opts.skips:
        pEnrichr.depends = pRNASeqDEG
        pEnrichr.args.pathview = False # don't do it for now
        pEnrichr.args.nthread = opts.nthread
        pEnrichr.args.inopts.cnames = True
        pEnrichr.args.genecol = opts.gscol - 1 if opts.gscol is not None else 0
        pEnrichr.args.libs = opts.enlibs
        pEnrichr.args.include = opts.enincl
        if opts.exdir:
            pEnrichr.config.export_dir = path.join(
                opts.exdir,
                '%d. Pathway enrichment analysis for DEGs' % dirindex
            )
            dirindex += 1

        pEnrichrUp = pEnrichr.copy() # pylint: disable=invalid-name
        pEnrichrUp.desc = 'Enrichment analysis for up-regulated genes'
        pEnrichrUp.depends = pRNASeqDEG
        pEnrichrUp.input = lambda ch: ch.outdir.expand(pattern='*.up.xls')
        if opts.exdir:
            pEnrichrUp.config.export_dir = path.join(
                opts.exdir,
                '%d. Pathway enrichment analysis for '
                'up-regulated DEGs' % dirindex
            )
            dirindex += 1

        pEnrichrDown = pEnrichr.copy() # pylint: disable=invalid-name
        pEnrichrDown.desc = 'Enrichment analysis for down-regulated genes'
        pEnrichrDown.depends = pRNASeqDEG
        pEnrichrDown.input = lambda ch: ch.outdir.expand(pattern='*.down.xls')
        if opts.exdir:
            pEnrichrDown.config.export_dir = path.join(
                opts.exdir,
                '%d. Pathway enrichment analysis for '
                'down-regulated DEGs' % dirindex
            )
            dirindex += 1

    ppl = PyPPL(ppldir=opts.ppldir, name="DEG_Analysis")
    ppl.start(starts)
    ppl.run(opts.runner)
    if opts.report:
        ppl.report(outfile=opts.report, standalone=False, title=opts.title)

if __name__ == "__main__":
    main()
