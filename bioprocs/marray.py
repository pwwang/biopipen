################################
#          processes           #
################################
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
pCeldir2Matrix.output              = "outfile:file:{{in.expdir | fn}}/{{in.expdir | fn}}.expr.txt, outdir:dir:{{in.expdir | fn}}"
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

################################
#         aggregations         #
################################
from pyppl import Aggr
from bioprocs.common import pPat2Dir, pFile2Proc
from bioprocs.resource import pTxt
from bioprocs.gsea import pExpmat2Gct, pSampleinfo2Cls, pGSEA, pEnrichr
"""
@name:
	aCelPat2Deg
@description:
	From celfils to degs with sample info file.
@input:
	`pattern`: The pattern to match the celfiles
	`sfile`  : The sample file
"""
aCelPat2Deg = Aggr(
	pPat2Dir,
	pFile2Proc,
	pCeldir2Matrix,
	pMarrayDeg,
	depends = False
)
# Dependences
aCelPat2Deg.starts                 = [aCelPat2Deg.pPat2Dir, aCelPat2Deg.pFile2Proc]
aCelPat2Deg.ends                   = [aCelPat2Deg.pMarrayDeg]
aCelPat2Deg.pCeldir2Matrix.depends = aCelPat2Deg.pPat2Dir
aCelPat2Deg.pMarrayDeg.depends     = aCelPat2Deg.pCeldir2Matrix, aCelPat2Deg.pFile2Proc
# Input
aCelPat2Deg.pMarrayDeg.input = lambda ch1, ch2: ch1.colAt(0).cbind(ch2)
# Args
aCelPat2Deg.pCeldir2Matrix.args.boxplot  = True
aCelPat2Deg.pCeldir2Matrix.args.heatmap  = True
aCelPat2Deg.pCeldir2Matrix.args.histplot = True
aCelPat2Deg.pMarrayDeg.args.maplot       = True
aCelPat2Deg.pMarrayDeg.args.heatmap      = True

"""
@name:
	aCelPat2DegGSEA
@description:
	From celfils to degs with sample info file and do GSEA.
@input:
	`pattern`: The pattern to match the celfiles
	`sfile`  : The sample file
	`gmtkey` : The gmtkey to gmt file to do the GSEA. See `bioprocs.resource.pTxt`
"""
aCelPat2DegGSEA = Aggr(
	pPat2Dir,
	pFile2Proc,
	pTxt,
	pCeldir2Matrix,
	pExpmat2Gct,
	pSampleinfo2Cls,
	pGSEA,
	pMarrayDeg,
	pEnrichr,
	depends = False
)
# Default input:
aCelPat2DegGSEA.pTxt.input = ['KEGG_2016_gmt']
# Dependences
aCelPat2DegGSEA.starts                  = aCelPat2DegGSEA.pPat2Dir, aCelPat2DegGSEA.pFile2Proc, aCelPat2DegGSEA.pTxt
aCelPat2DegGSEA.ends                    = aCelPat2DegGSEA.pGSEA, aCelPat2DegGSEA.pMarrayDeg, aCelPat2DegGSEA.pEnrichr
aCelPat2DegGSEA.pCeldir2Matrix.depends  = aCelPat2DegGSEA.pPat2Dir
aCelPat2DegGSEA.pExpmat2Gct.depends     = aCelPat2DegGSEA.pCeldir2Matrix
aCelPat2DegGSEA.pSampleinfo2Cls.depends = aCelPat2DegGSEA.pFile2Proc
aCelPat2DegGSEA.pGSEA.depends           = aCelPat2DegGSEA.pExpmat2Gct, aCelPat2DegGSEA.pSampleinfo2Cls, aCelPat2DegGSEA.pTxt
aCelPat2DegGSEA.pMarrayDeg.depends      = aCelPat2DegGSEA.pCeldir2Matrix, aCelPat2DegGSEA.pFile2Proc
aCelPat2DegGSEA.pEnrichr.depends        = aCelPat2DegGSEA.pMarrayDeg
# Input
aCelPat2DegGSEA.pMarrayDeg.input = lambda ch1, ch2: ch1.colAt(0).cbind(ch2)
# Args
aCelPat2DegGSEA.pCeldir2Matrix.args.boxplot  = True
aCelPat2DegGSEA.pCeldir2Matrix.args.heatmap  = True
aCelPat2DegGSEA.pCeldir2Matrix.args.histplot = True
aCelPat2DegGSEA.pMarrayDeg.args.maplot       = True
aCelPat2DegGSEA.pMarrayDeg.args.heatmap      = True
aCelPat2DegGSEA.pTxt.args.header             = False





