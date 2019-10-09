"""Microarray data analysis"""
import re
from os import path
from glob import glob
from pyppl import Proc, Box
#from .utils import plot, txt, dirnamePattern
from .rnaseq import pBatchEffect, pCoexp, pExprStats
from .utils import dirpat2name
from . import params, rimport
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pCELDir2Matrix():
	"""
	@name:
		pCELDir2Matrix
	@description:
		Convert CEL files to expression matrix
		File names will be used as sample names (colnames)
	@input:
		`indir:file`:  the directory containing the CEL files, could be gzipped
			- If you have files, then use `pFiles2Dir` first
		`sifile:File`: the sample infor file, extensions for samples are not necessary.
			- So that it can be also used by `pMArrayDEG`
	@output:
		`outfile:file`: the expression matrix file
		`outdir:dir`:   the directory containing expr file and plots
	@args:
		`pattern`  : The pattern to filter files. Default `'*'`
		`norm`     : The normalization method. Default: rma (mas5)
		`transfm`  : The extra tranformer for the expression values after the nomralization.
			- `Note the expression values have been done with log`
		`cdffile`  : The cdffile. Default: ''
		`fn2sample`: The transformer to transform file name
	"""
	pCELDir2Matrix                  = Proc(desc = 'Merge expression files to a matrix.')
	pCELDir2Matrix.input            = "indir:file, sifile:file"
	pCELDir2Matrix.output           = "outfile:file:{{i.indir, args.pattern | *dirpat2name}}.expr.txt"
	pCELDir2Matrix.lang             = params.Rscript.value
	pCELDir2Matrix.args.fn2sample   = 'function(fn) unlist(strsplit(fn, ".", fixed = T))[1]'
	pCELDir2Matrix.args.pattern     = '*'
	pCELDir2Matrix.args.norm        = 'rma' # mas5
	pCELDir2Matrix.args.transfm     = None # mas5
	pCELDir2Matrix.args.cdffile     = ''
	pCELDir2Matrix.envs.rimport     = rimport
	pCELDir2Matrix.envs.dirpat2name = dirpat2name
	pCELDir2Matrix.script           = "file:scripts/marray/pCELDir2Matrix.r"
	return pCELDir2Matrix

@procfactory
def _pMArrayDEG():
	"""
	@name:
		pMArrayDEG
	@description:
		Detect DEGs from microarray data.
	@input:
		`efile:file`: The expression matrix file
		`gfile:file`: The sample information file
	@output:
		`outfile:file`: The file with DEGs
		`outdir:dir`  : The directory containing files and plots
	@args:
		`tool`    : The tool to use. Default: `limma`
		`annofile`: The annotation file for the probes. Default: ``
			- If not provided, raw probe name will be used.
		`filter`: The filter of the data. Default: `[0, 0]`
			- 1st element: the `rowSums` of the expression matrix
			- 2nd element: how many samples(columns) have to reach the `rowSums` given by the 1st element
		`pval`  : The pval cutoff for DEGs. Default: `0.05`
			- You may specify which kind of measurement to use
			- `p:0.05`: using cutoff 0.05 for pvalue
			- `q:0.05`: using cutoff 0.05 for qvalue(fdr)
			- If no prefix is used, then it defaults to pvalue
		`plot`  : What kind of plots to generate.
			- `mdsplot` : The MDS plot
			- `volplot` : The volcano plot
			- `maplot`  : The MA plot
			- `heatmap` : The heatmap. (use `args.hmrows` to determine how many genes to plot)
		`ggs`   : The extra ggplot element for each plot (should be able to concatenate by `+`).
			- `maplot`  : for MA plot
			- `heatmap` : for heatmap. Default: `Box(theme = {'axis.text.y': 'r:element_blank()'})`
			- `volplot` : for volcano plot
		`devpars`: The parameters for plotting device. Default: `Box(res = 300, width = 2000, height = 2000)`
	@requires:
		`r-limma`
	"""
	return pMArrayDEG

@procfactory
def _pMArrayDEG():
	"""
	@name:
		pMArrayDEG
	"""
	pMArrayDEG        = Proc(desc = 'Detect DEGs from microarray data.')
	pMArrayDEG.input  = "efile:file, gfile:file"
	pMArrayDEG.output = [
		"outfile:file:{{i.efile | fn2}}-{{i.gfile | fn2}}-DEGs/{{i.efile | fn2}}-{{i.gfile | fn2}}.degs.txt",
		"outdir:dir:{{i.efile | fn2}}-{{i.gfile | fn2}}-DEGs"
	]
	pMArrayDEG.args.tool     = 'limma'
	pMArrayDEG.args.annofile = ''
	pMArrayDEG.args.filter   = [0, 0]
	pMArrayDEG.args.pval     = 0.05
	pMArrayDEG.args.hmrows   = 100
	pMArrayDEG.args.plot     = Box(
		mdsplot = True,
		volplot = Box(fccut = 2, pcut = 0.05),
		maplot  = False,
		heatmap = False
	)
	pMArrayDEG.args.ggs = Box(
		maplot  = Box(),
		heatmap = Box(theme = {'axis.text.y': 'r:element_blank()'}),
		volplot = Box(ylab = {0: '-log10(p-value)'})
	)
	pMArrayDEG.args.devpars = Box(res = 300, width = 2000, height = 2000)
	pMArrayDEG.envs.rimport = rimport
	pMArrayDEG.lang         = params.Rscript.value
	pMArrayDEG.script       = "file:scripts/marray/pMArrayDEG.r"
	return pMArrayDEG

