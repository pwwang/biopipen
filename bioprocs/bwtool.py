"""BigWig file operations using bwtool"""

from modkit import Modkit
from pyppl import Diot, Proc
from . import params, delefactory, procfactory
Modkit().delegate(delefactory())

@procfactory
def _pBWExtract():
	"""
	@description:
		Extract features from bigwig file using `bwtool extract`.
		See: https://github.com/CRG-Barcelona/bwtool/wiki/extract
	@input:
		infile: The input bigwig file.
		bedfile: The bed/gff file of regions to extract features from.
	@output:
		outfile: The output features from those regions
	@args:
		bwtool (str): Path to bwtool.
	"""
	return Diot(
		desc   = 'Extract features from bigwig file.',
		lang   = params.python.value,
		input  = 'infile:file, bedfile:file',
		output = 'outfile:file:{{i.infile | stem}}.bwtool-extracted.txt',
		args   = Diot(bwtool = params.bwtool.value)
	)

@procfactory
def _pBWSummary():
	"""
	@description:
		Provide some summary stats for each region in a bed file using `bwtool summary`
		See: https://github.com/CRG-Barcelona/bwtool/wiki/summary
		Only bed file for regions is supported in this process.
	@input:
		infile: The input bigwig file.
		bedfile: The bed/gff file of regions to extract features from.
	@output:
		outfile: The output file.
	@args:
		bwtool (str): The path to bwtool
		params (Diot): Other parameters for `bwtool summary`
	"""
	return Diot(
		desc   = 'provide some summary stats for each region in a bed file',
		lang   = params.python.value,
		input  = 'infile:file, bedfile:file',
		output = 'outfile:file:{{i.infile | stem}}.bwtool-summarized.txt',
		args   = Diot(
			bwtool = params.bwtool.value,
			params = Diot(header = True, keep_bed = True)
		)
	)
