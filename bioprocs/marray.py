from os import path
from glob import glob
from pyppl import Proc, Box
from .utils import plot, txt
from .rnaseq import pBatchEffect
from . import params

"""
@name:
	pCeldir2Matrix
@description:
	Convert CEL files to expression matrix
	File names will be used as sample names (colnames)
@input:
	`expdir:file`:  the directory containing the CEL files, could be gzipped
@output:
	`outfile:file`: the expression matrix file
	`outdir:dir`:   the directory containing expr file and plots
@args:
	`pattern` : The pattern to filter files. Default `'*'`
	`norm`    : The normalization method. Default: rma (mas5)
	`gfile`   : The group file. Default: ''
	`cdffile` : The cdffile. Default: ''
	`annofile`: The annotation file. Default: ''
	`exrows`  : Rows to be excluded, regular expression applied. Default: `[]`
	`boxplot` : Whether to plot a boxplot. Default: False
	`heatmap` : Whether to plot a heatmap. Default: False
	`histplot`: Whether to plot a histgram. Default: False
	`devpars` : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
	`boxplotggs`: The ggplot parameters for boxplot. Default: `['r:ylab("Expression")']`
		- See ggplot2 documentation.
	`heatmapggs`: The ggplot parameters for heatmap. Default: `['r:theme(axis.text.y = element_blank())']`
	`histplotggs`: The ggplot parameters for histgram. Default: `['r:labs(x = "Expression", y = "# Samples")']`	
"""
pCeldir2Matrix                     = Proc(desc = 'Merge expression files to a matrix.')
pCeldir2Matrix.input               = "expdir:file"
pCeldir2Matrix.output              = [
	"outfile:file:{{in.expdir, args.pattern | fsDirname}}/{{in.expdir, args.pattern | fsDirname}}.expr.txt", 
	"outdir:dir:{{in.expdir, args.pattern | fsDirname}}"
]
pCeldir2Matrix.lang                = params.Rscript.value
pCeldir2Matrix.args.pattern        = '*'
pCeldir2Matrix.args.norm           = 'rma' # mas5
pCeldir2Matrix.args.gfile          = ''
pCeldir2Matrix.args.cdffile        = ''
pCeldir2Matrix.args.annofile       = ''
pCeldir2Matrix.args.boxplot        = False
pCeldir2Matrix.args.heatmap        = False
pCeldir2Matrix.args.heatmapn       = 500
pCeldir2Matrix.args.histplot       = False
pCeldir2Matrix.args.devpars        = Box({'res': 300, 'width': 2000, 'height': 2000})
pCeldir2Matrix.args.boxplotggs     = ['r:ylab("Log2 Intensity")']
pCeldir2Matrix.args.heatmapggs     = ['r:theme(axis.text.y = element_blank())']
pCeldir2Matrix.args.histplotggs    = ['r:labs(x = "Log2 Intensity", y = "Density")']
pCeldir2Matrix.tplenvs.plotBoxplot = plot.boxplot.r
pCeldir2Matrix.tplenvs.plotHeatmap = plot.heatmap.r
pCeldir2Matrix.tplenvs.plotHist    = plot.hist.r
pCeldir2Matrix.tplenvs.fsDirname   = lambda dir, pat: path.splitext(path.basename(glob(path.join(dir, pat))[0]))[0] + '_etc'
pCeldir2Matrix.script              = "file:scripts/marray/pCeldir2Matrix.r"



pMarrayDeg        = Proc(desc = 'Detect DEGs by microarray data.')
pMarrayDeg.input  = "efile:file, gfile:file"
pMarrayDeg.output = [
	"outfile:file:{{in.efile | fn | fn}}-{{in.gfile | fn | fn}}-DEGs/{{in.efile | fn | fn}}-{{in.gfile | fn | fn}}.degs.txt", 
	"outdir:dir:{{in.efile | fn | fn}}-{{in.gfile | fn | fn}}-DEGs"
]
pMarrayDeg.args.tool             = 'limma'
pMarrayDeg.args.filter           = '1,2'
pMarrayDeg.args.pval             = 0.05
pMarrayDeg.args.mdsplot          = True
pMarrayDeg.args.volplot          = True
pMarrayDeg.args.maplot           = False
pMarrayDeg.args.heatmap          = False
pMarrayDeg.args.heatmapn         = 100
pMarrayDeg.args.heatmapggs       = ['r:theme(axis.text.y = element_blank())']
pMarrayDeg.args.maplotggs        = []
pMarrayDeg.args.volplotggs       = []
pMarrayDeg.args.devpars          = Box({'res': 300, 'width': 2000, 'height': 2000})
pMarrayDeg.tplenvs.plotHeatmap   = plot.heatmap.r
pMarrayDeg.tplenvs.plotMAplot    = plot.maplot.r
pMarrayDeg.tplenvs.plotVolplot   = plot.volplot.r
pMarrayDeg.tplenvs.txtSampleinfo = txt.sampleinfo.r
pMarrayDeg.lang                  = params.Rscript.value
pMarrayDeg.script                = "file:scripts/marray/pMarrayDeg.r"
