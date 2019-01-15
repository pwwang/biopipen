"""
Generate plot using given data
"""
from pyppl import Proc, Box
#from .utils import plot
from . import params, rimport

"""
@name:
	pPlot
@description:
	Use ggplot2 to generate plots
@input:
	`infile:file`: The input data file
@output:
	`outfile:file`: The output file
@args:
	`cnames` : Whether the input file has colnames. Default: True
	`rnames` : Whether the input file has rownames. Default: False
	`aes`    : The default aes. Default: {'x':1, 'y':2} (corresponding to colnames)
	`helper` : Some helper codes to generate `params` and `ggs`
	`devpars`: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`
	`ggs`    : The extra ggplot elements.
"""
pPlot = Proc(desc = 'Generate plot using ggplot2')
pPlot.input  = 'infile:file'
pPlot.output = 'outfile:file:{{i.infile | fn}}.png'
pPlot.args.cnames  = True
pPlot.args.rnames  = False
pPlot.args.aes     = Box(x = 1, y = 2) # only allow x, y
pPlot.args.helper  = ''
pPlot.args.devpars = Box(res = 300, height = 2000, width = 2000)
pPlot.args.params  = Box()
pPlot.args.ggs     = Box()
pPlot.envs.rimport = rimport
pPlot.lang   = params.Rscript.value
pPlot.script = 'file:scripts/plot/pPlot.r'

"""
@name:
	pScatter
@description:
	Use ggplot2 geom_point to generate plots
@infile:
	`infile:file`: The input data file
@outfile:
	`outfile:file`: The output file
@args:
	`cnames` : Whether the input file has colnames. Default: True
	`rnames` : Whether the input file has rownames. Default: False
	`x`      : The x aes. Default: 1 (corresponding to colnames)
	`y`      : The y aes. Default: 2 (corresponding to colnames)
	`helper` : Some helper codes to generate `params` and `ggs`
	`devpars`: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`
	`params` : The extra params for `geom_point`
	`ggs`    : The extra ggplot elements.
"""
pScatter              = Proc(desc = 'Generate scatter plot.')
pScatter.input        = 'infile:file'
pScatter.output       = 'outfile:file:{{i.infile | fn}}.scatter.png'
pScatter.args.cnames  = True
pScatter.args.rnames  = False
pScatter.args.x       = 1
pScatter.args.y       = 2
pScatter.args.helper  = ''
pScatter.args.devpars = Box(res = 300, height = 2000, width = 2000)
pScatter.args.params  = Box()
pScatter.args.ggs     = Box()
pScatter.envs.rimport = rimport
pScatter.lang         = params.Rscript.value
pScatter.script       = 'file:scripts/plot/pScatter.r'

"""
@name:
	pPoints
@description:
	Alias for pScatter
"""
pPoints = pScatter.copy()

"""
@name:
	pHisto
@description:
	Use ggplot2 geom_histogram to generate histograms
@infile:
	`infile:file`: The input data file
@outfile:
	`outfile:file`: The output file
@args:
	`cnames` : Whether the input file has colnames. Default: True
	`rnames` : Whether the input file has rownames. Default: False
	`x`      : The x aes. Default: 1 (corresponding to colnames)
	`helper` : Some helper codes to generate `params` and `ggs`
	`devpars`: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`
	`params` : The extra params for `geom_point`
	`ggs`    : The extra ggplot elements.
"""
pHisto              = Proc(desc = 'Generate histogram.')
pHisto.input        = 'infile:file'
pHisto.output       = 'outfile:file:{{i.infile | fn}}.histo.png'
pHisto.args.cnames  = True
pHisto.args.rnames  = False
pHisto.args.x       = 1
pHisto.args.helper  = ''
pHisto.args.devpars = Box(res = 300, height = 2000, width = 2000)
pHisto.args.params  = Box()
pHisto.args.ggs     = Box()
pHisto.envs.rimport = rimport
pHisto.lang         = params.Rscript.value
pHisto.script       = 'file:scripts/plot/pHisto.r'

"""
@name:
	pFreqpoly
@description:
	Use ggplot2 geom_freqpoly to generate frequency polygon plot.
@infile:
	`infile:file`: The input data file
@outfile:
	`outfile:file`: The output file
@args:
	`cnames` : Whether the input file has colnames. Default: True
	`rnames` : Whether the input file has rownames. Default: False
	`x`      : The x aes. Default: 1 (corresponding to colnames)
	`helper` : Some helper codes to generate `params` and `ggs`
	`devpars`: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`
	`params` : The extra params for `geom_point`
	`ggs`    : The extra ggplot elements.
"""
pFreqpoly = Proc(desc = 'Generate frequency polygon plot.')
pFreqpoly.input  = 'infile:file'
pFreqpoly.output = 'outfile:file:{{i.infile | fn}}.freqpoly.png'
pFreqpoly.args.cnames  = True
pFreqpoly.args.rnames  = False
pFreqpoly.args.x       = 1
pFreqpoly.args.helper  = ''
pFreqpoly.args.devpars = Box(res = 300, height = 2000, width = 2000)
pFreqpoly.args.params  = Box()
pFreqpoly.args.ggs     = Box()
pFreqpoly.envs.rimport = rimport
pFreqpoly.lang   = params.Rscript.value
pFreqpoly.script = 'file:scripts/plot/pFreqpoly.r'


"""
@name:
	pBoxplot
@description:
	Generate box plot
@input:
	`datafile:file`: The data file
@output:
	`outpng:file`: The output figure
@args:
	`inopts` :   Input options to read the input file
		- `cnames` :   Whether the input file has header. Default: `True`
		- `rnames` :   Whether the input file has row names. Default: `False`
		- `delimit`:   The seperator. Defualt: `\\t`
	`x`      :   The `ind` (index) column. Only for `args.stacked = True`. Default: `2`
	`y`      :   The `values` column. Only for `args.stacked = True`. Default: `1`
	`helper` :   Some raw codes to help to construct the matrix and arguments.
	`stacked`:   Whether the input file is stacked
		- Stacked file looks like:
		  ```
		  values	ind
		  1.1	col1
		  1.2	col1
		  ...
		  .8	col2
		  .9	col2
		  ...
		  3.2	col3
		  ...
		  ```
		- Unstacked file looks like:
		  ```
		  col1	col2	col3
		  1.1	.8	3.2
		  1.2	.9	2.2
		  ```
	`params`:    Other parameters for `boxplot`, default: `""`
	`ggs`   :    Extra ggplot2 statements
"""
pBoxplot             = Proc(desc = 'Generate boxplot plot.')
pBoxplot.input       = 'infile:file'
pBoxplot.output      = 'outfile:file:{{i.infile | fn}}.boxplot.png'
pBoxplot.args.inopts = Box(
	cnames  = True,
	rnames  = False,
	delimit = '\t'
)
pBoxplot.args.x       = 2
pBoxplot.args.y       = 1
pBoxplot.args.helper  = ''
pBoxplot.args.stacked = False
pBoxplot.args.devpars = Box(res = 300, height = 2000, width = 2000)
pBoxplot.args.params  = Box()
pBoxplot.args.ggs     = Box()
pBoxplot.envs.rimport = rimport
pBoxplot.lang         = params.Rscript.value
pBoxplot.script       = 'file:scripts/plot/pBoxplot.r'

"""
@name:
	pHeatmap
@description:
	Plot heatmaps.
@input:
	`infile:file`: The input matrix file
@output:
	`outfile:file`: The heatmap
@args:
	`ggs`: The ggplot items for heatmap
	`devpars`: The parameters for device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`
	`dendro`: The parameters for control of the dendrogram. Default: `{'dendro': True}`
		- `dendro`: `True`: plot dendros for both rows and cols; `col`: only plot dendro for cols; `row`: only plot dendro for rows
		- `rows`: The rownames to subset the rows and control the order of rows. Must a list. Only works when not plotting dendro for rows.
		- `cols`: The colnames to subset the cols and control the order of cols. Must a list. Only works when not plotting dendro for cols.
	`header`: The input file has header? Default: True
	`rownames`: The input file has rownames? Default: 1
	`rows`: Row selector
		- `all`: All rows
		- `top:N`: Top N rows (original data ordered in descending order). N defaults to 100
		- `bottom:N`: Bottom N rows. N defaults to 100
		- `both:N`: Top N rows and bottom N rows. N defaults to 50
		- `random:N`: Random N rows. N defaults to 50
		- `random-both:N`: Random N rows from top part and N rows from bottom part. N defaults to 50
	`cols`: Col selector (see `rows`).
"""
pHeatmap              = Proc(desc  = 'Plot heatmaps.')
pHeatmap.input        = "infile:file"
pHeatmap.output       = "outfile:file:{{i.infile | fn}}.heatmap.png"
pHeatmap.args.ggs     = Box()
pHeatmap.args.devpars = Box(res = 300, height = 2000, width = 2000)
pHeatmap.args.params  = Box(dendro = True)
pHeatmap.args.inopts  = Box(rnames = True, cnames = True)
pHeatmap.args.helper  = ''
pHeatmap.envs.rimport = rimport
pHeatmap.lang         = params.Rscript.value
pHeatmap.script       = "file:scripts/plot/pHeatmap.r"

"""
@name:
	pScatterCompare
@description:
	Plot scatter plot to compare values of first 2 columns of input data
@input:
	`infile:file`: The input file containing a matrix with at least 2 columns
		- Other columns are groups used to group the scatter points
		- Data must be normalized to [0, 1]
@output:
	`outfile:file`: The output plot
@args:
	`ggs`: Extra expressions for ggplot. Note if geom_point is included, original geom_point will be ignored.
	`devpars`: The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`
	`rownames`: Whether the input file has row names. Default: True
	`regr`: Whether draw the regression line. Default: False
	`corr`: The method to calculate the correlation. Default: `pearson`
		- Could be: `pearson`, `spearman` or `kendall`
		- If it's neither of the three, no correlations will show.
"""
pScatterCompare              = Proc(desc = 'Plot scatter compare plots.')
pScatterCompare.input        = "infile:file"
pScatterCompare.output       = "outfile:file:{{i.infile | fn}}.scattercomp.png"
pScatterCompare.args.ggs     = Box()
pScatterCompare.args.params  = Box()
pScatterCompare.args.devpars = Box(res = 300, height = 2000, width = 2000)
pScatterCompare.args.x       = 1
pScatterCompare.args.y       = 2
pScatterCompare.args.rnames  = True
pScatterCompare.args.cnames  = True
pScatterCompare.args.helper  = ''
pScatterCompare.envs.rimport = rimport
pScatterCompare.lang         = params.Rscript.value
pScatterCompare.script       = "file:scripts/plot/pScatterCompare.r"

"""
@name:
	pROC
@description:
	Generate ROC curves and output AUC.
@input:
	`infile:file`: The input matrix file.
		- Col1: rownames if args.rnames is True else label (0, 1 class)
		- Col2: prediction values from model1
		- ...
@output:
	`outdir:dir`: The output directory
@args:
	`inopts`: The options for input file. Default: `Box(rnames = True, cnames = True)`
	`params`: The parameters for `plot.roc` from `utils/plot.r`
	`ggs`   : Additaional ggplot terms. Default: 
		```python
		Box({
		    'style_roc': {},
		    # show legend at bottom right corner
		    'theme#auc': {'legend.position': [1, 0], 'legend.justification': [1, 0]} 
		})
		```
	`devpars`: The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`
"""
pROC              = Proc(desc = 'Generate ROC curves.')
pROC.input        = 'infile:file'
pROC.output       = [
	'outfile:file:{{i.infile | fn}}.roc/{{i.infile | fn}}.auc.txt', 
	'outdir:dir:{{i.infile | fn}}.roc'
]
pROC.args.inopts  = Box(rnames = True, cnames = True)
pROC.args.params  = Box(labels = False, showAUC = True, combine = True)
pROC.args.ggs     = Box({
	'style_roc': {},
	# show legend at bottom right corner
	'theme#auc': {'legend.position': [1, 0], 'legend.justification': [1, 0]} 
})
pROC.args.devpars = Box(res = 300, height = 2000, width = 2000)
pROC.envs.rimport = rimport
pROC.lang         = params.Rscript.value
pROC.script       = "file:scripts/plot/pROC.r"

"""
@name:
	pVenn
@description:
	Venn/UpsetR plots.
@input:
	`infile:file`: The input matrix
		- format:
		```
			category1	category2	category3
		[e1]	0	1	1
		[e2]	0	0	1
		...
		[eN]	1	0	0
		```
		rownames are not necessary but colnames are.
	`metafile:file`: The metadata file for each category for upset plot.
		- format:
		```
			col1	col2	...	colN
		category1	x1	y1	...	z1
		category2	x2	y2	...	z2
		...	...
		categoryN	xN	yN	...	zN
		```
@output:
	`outfile:file`: The plot
@args:
	`tool`    : Which tools to use. Default: auto (venn, upsetr, auto(n<=3: venn, otherwise upsetr))
	`rnames`  : Whether input file has rownames. Default: False
	`params`  : Other params for `venn.diagram` or `upset`. Default: {}
	`devpars` : The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`
@requires:
	[`r-VennDiagram`](https://www.rdocumentation.org/packages/VennDiagram)
	[`r-UpSetR`](https://www.rdocumentation.org/packages/UpSetR)
"""
pVenn              = Proc(desc = 'Venn plots.')
pVenn.input        = "infile:file, metafile:file"
pVenn.output       = "outfile:file:{{i.infile | fn}}.venn.png"
pVenn.args.tool    = 'auto' # upsetr or auto: < = 3 venn, else upsetr
pVenn.args.rnames  = False
pVenn.args.params  = Box()
pVenn.args.devpars = Box(res = 300, height = 2000, width = 2000)
pVenn.envs.rimport = rimport
pVenn.lang         = params.Rscript.value
pVenn.script       = "file:scripts/plot/pVenn.r"

"""
@name:
	pPie
@description:
	Plot piechart
@input:
	`infile:file`: The input file. Could be either:
		- Direct numbers of each category.
		```
		Group1	Group2
		50	50
		```
		- Presence of each items in the category.
		```
			Group1	Group2
		Item1	1	0
		Item2	0	1
		...
		```
@output:
	`outfile:file`: the output plot
@args:
	`rnames` : Whether the input file has row names. Default: `False`
	`ggs`    : Extra expressions for ggplot.
	`devpars`: The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`
"""
pPie              = Proc(desc = 'Pie chart.')
pPie.input        = "infile:file"
pPie.output       = "outfile:file:{{i.infile | fn}}.pie.png"
pPie.args.rnames  = False
pPie.args.devpars = Box(res = 300, height = 2000, width = 2000)
pPie.args.ggs     = Box(
	theme = {
		'axis.title.x'    : 'r:element_blank()',
		'axis.title.y'    : 'r:element_blank()',
		'panel.border'    : 'r:element_blank()',
		'panel.grid'      : 'r:element_blank()',
		'axis.ticks'      : 'r:element_blank()',
		'panel.background': 'r:element_blank()',
		'axis.text.x'     : 'r:element_blank()',
	}
)
pPie.envs.rimport = rimport
pPie.lang         = params.Rscript.value
pPie.script       = "file:scripts/plot/pPie.r"

pManhattan              = Proc(desc = 'Manhattan plot.')
pManhattan.input        = 'infile:file, hifile:file'
pManhattan.output       = 'outfile:file:{{i.infile | fn}}.manht.png'
pManhattan.args.inopts  = Box(cnames = False, rnames = False)
pManhattan.args.devpars = Box(res = 300, height = 2000, width = 2000)
pManhattan.args.ggs     = Box()
pManhattan.envs.rimport = rimport
pManhattan.lang         = params.Rscript.value
pManhattan.script       = "file:scripts/plot/pManhattan.r"

