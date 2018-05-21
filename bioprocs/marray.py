import re
from os import path
from glob import glob
from pyppl import Proc, Box
#from .utils import plot, txt, dirnamePattern
from .rnaseq import pBatchEffect, pCoexp
from .utils import dirpat2name
from . import params, rimport

"""
@name:
	pCELdir2Matrix
@description:
	Convert CEL files to expression matrix
	File names will be used as sample names (colnames)
@input:
	`indir:file`:  the directory containing the CEL files, could be gzipped
		- If you have files, then use `pFiles2Dir` first
@output:
	`outfile:file`: the expression matrix file
	`outdir:dir`:   the directory containing expr file and plots
@args:
	`pattern`  : The pattern to filter files. Default `'*'`
	`norm`     : The normalization method. Default: rma (mas5)
	`gfile`    : The group file. Default: ''
	`cdffile`  : The cdffile. Default: ''
	`annofile` : The annotation file. Default: ''
	`hmrows`   : How many rows to be used to plot heatmap
	`plot`: Whether to plot
		- `boxplot`   : Whether to plot a boxplot. Default: False
		- `heatmap`   : Whether to plot a heatmap. Default: False
		- `histogram` : Whether to plot a histgram. Default: False
	`devpars`    : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
	`ggs`: The ggplot parameters
		- `boxplot`  : The ggplot parameters for boxplot. Default: `Box(ylab = {0: "Log2 Intensity"})`
		- `heatmap`  : The ggplot parameters for heatmap. Default: `Box(theme = {'axis.text.y': 'r:element_blank()'})`
		- `histogram`: The ggplot parameters for histgram. Default: `Box(labs = {'x': "Log2 Intensity", "y": "Density"})`
"""
pCELdir2Matrix               = Proc(desc = 'Merge expression files to a matrix.')
pCELdir2Matrix.input         = "indir:file"
pCELdir2Matrix.output        = [
	"outfile:file:{{in.indir, args.pattern | dirpat2name}}.dir/{{in.indir, args.pattern | dirpat2name}}.expr.txt",
	"outdir:dir:{{in.indir, args.pattern | dirpat2name}}.dir"
]
pCELdir2Matrix.lang           = params.Rscript.value
pCELdir2Matrix.args.fn2sample = 'function(fn) fn'
pCELdir2Matrix.args.pattern   = '*'
pCELdir2Matrix.args.norm      = 'rma' # mas5
pCELdir2Matrix.args.gfile     = ''
pCELdir2Matrix.args.cdffile   = ''
pCELdir2Matrix.args.annofile  = ''
pCELdir2Matrix.args.hmrows    = 500
pCELdir2Matrix.args.plot      = Box(boxplot = False, heatmap = False, histogram = False)
pCELdir2Matrix.args.ggs       = Box(
	boxplot   = Box(ylab  = {0: "Log2 Intensity"}),
	heatmap   = Box(theme = {'axis.text.y': 'r:element_blank()'}),
	histogram = Box(labs  = {'x': "Log2 Intensity", "y": "Density"})
)
pCELdir2Matrix.args.devpars     = Box(res = 300, width = 2000, height = 2000)
pCELdir2Matrix.envs.rimport     = rimport
pCELdir2Matrix.envs.dirpat2name = dirpat2name
pCELdir2Matrix.script           = "file:scripts/marray/pCELdir2Matrix.r"



pMArrayDEG        = Proc(desc = 'Detect DEGs by microarray data.')
pMArrayDEG.input  = "efile:file, gfile:file"
pMArrayDEG.output = [
	"outfile:file:{{in.efile | fn2}}-{{in.gfile | fn2}}-DEGs/{{in.efile | fn2}}-{{in.gfile | fn2}}.degs.txt",
	"outdir:dir:{{in.efile | fn2}}-{{in.gfile | fn2}}-DEGs"
]
pMArrayDEG.args.tool   = 'limma'
pMArrayDEG.args.filter = '1,2'
pMArrayDEG.args.pval   = 0.05
pMArrayDEG.args.hmrows = 100
pMArrayDEG.args.plot   = Box(
	mdsplot = True,
	volplot = True,
	maplot  = False,
	heatmap = False
)
pMArrayDEG.args.ggs = Box(
	maplot  = Box(),
	heatmap = Box(theme = {'axis.text.y': 'r:element_blank()'}),
	volplot = Box()
)
pMArrayDEG.args.devpars = Box(res = 300, width = 2000, height = 2000)
pMArrayDEG.envs.rimport = rimport
pMArrayDEG.lang         = params.Rscript.value
pMArrayDEG.script       = "file:scripts/marray/pMArrayDEG.r"
