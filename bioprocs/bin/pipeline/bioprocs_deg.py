#!/usr/bin/env python
"""Run DEG analysis for expression data"""

# Core information needed:
# 1. expression matrix, genes as rows and samples as columns, must be count data.
# 2. sample information file, specifying Samples, Patients (if paired analysis) and Group (comparison)

from os import path

from pyppl import PyPPL, Box
from bioprocs import params

params._prefix = '-'
params._desc = [
	"Things this script can do:",
	"1. Statistics/Characteristics of the expression matrix",
	"2. DEG Calling",
	"3. Real GSEA analysis",
	"4. Enrichment (depends on step 2)"]
params.exprmat.required = True
params.exprmat.desc     = 'The expression matrix.'
params.saminfo.required = True
params.saminfo.desc     = 'The sample information file.'
params.ppldir           = './workdir'
params.ppldir.desc      = 'The pipeline directory.'
params.exdir.required   = True
params.exdir.desc       = 'Where to export the result files.'
params.caller           = 'deseq2'
params.caller.desc      = 'The DEG caller'
params.runner           = 'local'
params.runner.desc      = 'The runner to run the processes.'
params.mapping          = ''
params.mapping.desc     = 'Probe to gene mapping'
params.nthread          = 1
params.nthread.desc     = '# threads to use.'
params.enlibs           = ['KEGG_2016']
params.enlibs.desc      = 'The enrichment analysis libraries'
params.cutoff           = 0.05
params.cutoff.desc      = 'The cutoff for DEGs'
params.cutoff.callback  = lambda opt: setattr(opt, 'value', {"by": "q", "value": opt.value} if not isinstance(opt.value, dict) else opt.value)
params.filter           = ""
params.filter.desc      = "A filter for the expression matrix"
params.filter.callback  = lambda opt: setattr(opt, 'value', opt.value or None)
params.steps            = ['stats', 'call', 'gsea', 'enrich']
params.steps.desc       = [
	'The steps to do with this script',
	'Use "-step" to clear defaults, and "-step <STEP>" to select step',
	'For example: "-step -step stats" will only run "stats",',
	'Only "-step stats" will not clear the default steps.'
]
params.steps.callback   = lambda opt: all(v in ['stats', 'call', 'gsea', 'enrich'] for v in opt.value) or "-steps should be a subset of ['stats', 'call', 'gsea', 'enrich']"

def main():
	opts = params._parse(dict_wrapper = Box)

	from bioprocs.gsea import pEnrichr
	from bioprocs.rnaseq import pExprStats, pRNASeqDEG

	# processes
	pExprStats.input = (opts.exprmat, opts.saminfo)
	pExprStats.args.tsform = 'function(x) log2(x+1)'
	pExprStats.args.opts.histogram.bins = 100

	pRNASeqDEG.args.tool    = opts.caller
	pRNASeqDEG.exdir        = path.join(opts.exdir, '2.DEGs')
	pRNASeqDEG.args.cutoff  = opts.cutoff
	pRNASeqDEG.args.mapfile = opts.mapping

	pEnrichr.args.genecol = 1 if opts.mapping else 0
	pEnrichr.args.libs = opts.enlibs
	pEnrichr.exdir = path.join(opts.exdir, '3.Enrichment', 'AllDEGs')

	if opts.filter:
		pExprStatsFiltered = pExprStats.copy()
		pExprStats.exdir = path.join(opts.exdir, '0.Unfiltered.expression.stats')
		pExprStatsFiltered.exdir = path.join(opts.exdir, '1.Filtered.expression.stats')
		pExprStatsFiltered.args.filter = opts.filter
	else:
		pExprStats.exdir = path.join(opts.exdir, '1.Expression.stats')

	# starts/depends
	starts = []
	if 'stats' in opts.steps:
		starts.append(pExprStats)
		if opts.filter:
			starts.append(pExprStatsFiltered)
	if 'call' in opts.steps:
		if 'stats' in opts.steps and opts.filter:
			pRNASeqDEG.depends = pExprStatsFiltered
			pRNASeqDEG.input = lambda ch: ch.expand(pattern = path.basename(opts.exprmat)).cbind(opts.saminfo)
		else:
			pRNASeqDEG.input = (opts.exprmat, opts.saminfo)
			starts.append(pRNASeqDEG)
	if 'gsea' in opts.steps:
		# to be implemented
		pass
	if 'enrich' in opts.steps:
		if 'call' not in opts.steps:
			raise ValueError('"enrich" step requires "call" step.')
		pEnrichr.depends             = pRNASeqDEG
		pEnrichr.args.pathview.fccol = int(opts.caller == 'deseq2') + int(bool(opts.mapping)) + 2
		pEnrichr.args.nthread        = opts.nthread
		pEnrichr.args.genecol        = int(bool(opts.mapping))
		pEnrichrUp                   = pEnrichr.copy()
		pEnrichrUp.exdir             = path.join(opts.exdir, '3.Enrichment', 'UpDEGs')
		pEnrichrUp.depends           = pRNASeqDEG
		pEnrichrUp.input             = lambda ch: ch.outdir.expand(pattern = '*.up.txt')
		pEnrichrDown                 = pEnrichr.copy()
		pEnrichrDown.exdir           = path.join(opts.exdir, '3.Enrichment', 'DownDEGs')
		pEnrichrDown.depends         = pRNASeqDEG
		pEnrichrDown.input           = lambda ch: ch.outdir.expand(pattern = '*.down.txt')

	config = {
		'default': {'ppldir': opts.ppldir, 'acache': True}
	}
	PyPPL(config).start(starts).run(opts.runner)

if __name__ == "__main__":
	main()


