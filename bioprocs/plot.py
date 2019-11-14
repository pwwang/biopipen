"""Generate plots using given data"""

from pyppl import Proc, Box
from . import params, rimport
from .utils import fs2name
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pPlot():
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
	return pPlot

@procfactory
def _pScatter(alias = 'pPoints'):
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
	return pScatter

@procfactory
def _pHisto():
	"""
	@name:
		pHisto
	@description:
		Use ggplot2 geom_histogram to generate histograms
	@input:
		`infile:file`: The input data file
	@output:
		`outfile:file`: The output file
	@args:
		`inopts` : Options to read the input file. Default: `Box(cnames = True, rnames = False)`
		`x`      : The x aes. Default: 1 (corresponding to colnames)
		`helper` : Some helper codes to generate `params` and `ggs`
		`devpars`: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`
		`params` : The extra params for `geom_histogram`
		`ggs`    : The extra ggplot elements.
	"""
	pHisto              = Proc(desc = 'Generate histogram.')
	pHisto.input        = 'infile:file'
	pHisto.output       = 'outfile:file:{{i.infile | fn}}.histo.png'
	pHisto.args.inopts  = Box(cnames = True, rnames = False)
	pHisto.args.x       = 1
	pHisto.args.helper  = ''
	pHisto.args.devpars = Box(res = 300, height = 2000, width = 2000)
	pHisto.args.params  = Box()
	pHisto.args.ggs     = Box()
	pHisto.envs.rimport = rimport
	pHisto.lang         = params.Rscript.value
	pHisto.script       = 'file:scripts/plot/pHisto.r'
	return pHisto

@procfactory
def _pDensity():
	"""
	@name:
		pDensity
	@description:
		Use ggplot2 geom_density to generate density plot
	@input:
		`infile:file`: The input data file
	@output:
		`outfile:file`: The output file
	@args:
		`inopts` : Options to read the input file. Default: `Box(cnames = True, rnames = False)`
		`x`      : The x aes. Default: 1 (corresponding to colnames)
		`helper` : Some helper codes to generate `params` and `ggs`
		`devpars`: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`
		`params` : The extra params for `geom_density`
		`ggs`    : The extra ggplot elements.
	"""
	pDensity              = Proc(desc = 'Generate density plot.')
	pDensity.input        = 'infile:file'
	pDensity.output       = 'outfile:file:{{i.infile | fn}}.density.png'
	pDensity.args.inopts  = Box(cnames = True, rnames = False)
	pDensity.args.x       = 1
	pDensity.args.helper  = ''
	pDensity.args.devpars = Box(res = 300, height = 2000, width = 2000)
	pDensity.args.params  = Box()
	pDensity.args.ggs     = Box()
	pDensity.envs.rimport = rimport
	pDensity.lang         = params.Rscript.value
	pDensity.script       = 'file:scripts/plot/pDensity.r'
	return pDensity

@procfactory
def _pFreqpoly():
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
	return pFreqpoly


@procfactory
def _pBoxplot():
	"""
	@name:
		pBoxplot
	@description:
		Generate box plot
	@input:
		`infile:file`: The data file
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
		`params`:    Other parameters for `geom_boxplot`, default: `Box()`
		`ggs`   :    Extra ggplot2 statements
	"""
	pBoxplot              = Proc(desc = 'Generate boxplot plot.')
	pBoxplot.input        = 'infile:file'
	pBoxplot.output       = 'outfile:file:{{i.infile | fn}}.boxplot.png'
	pBoxplot.args.inopts  = Box(cnames = True, rnames = False, delimit = '\t')
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
	return pBoxplot

@procfactory
def _pBar(alias = 'pCol'):
	"""
	@name:
		pBar
	@description:
		Generate bar/col plot
	@input:
		`infile:file`: The data file
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
			- see `pBoxplot.args.stacked`
		`params`:    Other parameters for `geom_bar`, default: `Box()`
		`ggs`   :    Extra ggplot2 statements
	"""
	pBar             = Proc(desc = 'Generate bar/col plot.')
	pBar.input       = 'infile:file'
	pBar.output      = 'outfile:file:{{i.infile | fn}}.bar.png'
	pBar.args.inopts = Box(
		cnames  = True,
		rnames  = False,
		delimit = '\t'
	)
	pBar.args.x       = 2
	pBar.args.y       = 1
	pBar.args.helper  = ''
	pBar.args.stacked = False
	pBar.args.devpars = Box(res = 300, height = 2000, width = 2000)
	pBar.args.params  = Box()
	pBar.args.ggs     = Box()
	pBar.envs.rimport = rimport
	pBar.lang         = params.Rscript.value
	pBar.script       = 'file:scripts/plot/pBar.r'
	return pBar

@procfactory
def _pHeatmap():
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
	return pHeatmap

@procfactory
def _pHeatmap2():
	"""
	@name:
		pHeatmap2
	@description:
		Plot heatmaps using R package ComplexHeatmap. Example:
		```
		bioprocs plot.pHeatmap2
		-i.infile MMPT.txt
		-i.annofiles:l:o PatientAnno.txt
		-args.params.row_names_gp 'r:fontsize5'
		-args.params.column_names_gp 'r:fontsize5'
		-args.params.clustering_distance_rows pearson
		-args.params.clustering_distance_columns pearson
		-args.params.show_row_names false
		-args.params.row_split 3
		-args.devpars.width 5000
		-args.devpars.height 5000
		-args.draw.merge_legends
		-args.params.heatmap_legend_param.title AUC
		-args.params.row_dend_reorder
		-args.params.column_dend_reorder
		-args.params.top_annotation 'r:HeatmapAnnotation(Mutation = as.matrix(annos[,(length(groups)+1)  :ncol(annos)]), Group = as.matrix(annos[,groups]), col = list(Mutation = c(`0`="grey",   `1`="lightgreen", `2`="green", `3`="darkgreen")), annotation_name_gp = fontsize8, show_legend =   c(Group=F))'
		-args.params.right_annotation 'r:rowAnnotation(AUC = anno_boxplot(as.matrix(data), outline = F))  '
		-args.helper 'fontsize8 = gpar(fontsize = 12); fontsize5 = gpar(fontsize = 8); groups = c  ("Group1", "Group2", "Group3")'
		-args.seed 8525
		```
	@input:
		`infile:file`: The input data file for the main heatmap.
		`annofiles:files`: The annotation files.
			- For now, they should share the same `args.anopts`
	@output:
		`outfile:file`: The plot.
		`outdir:dir`: The output directory including the output plot and other information.
	@args:
		`devpars` : The parameters for device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`
		`draw`    : The parameters for `ComplexHeatmap::draw`
		`params`  : Other parameters for `ComplexHeatmap::Heatmap`
		`anopts`  : The options to read the annotation files.
		`inopts`  : The options to read the input files.
		`seed`    : The seed. Default: `None`
		`helper`  : The raw R codes to help defining some R variables or functions.
		`saveinfo`: Save other information include the heatmap object and the clusters. Default: `True`
	@requires:
		`R-ComplexHeatmap` (tested on c269eb425bf1b2d1713b9d5e68bf9f08bd8c7acb)
	"""
	pHeatmap2              = Proc(desc  = 'Plot heatmaps using R package ComplexHeatmap.')
	pHeatmap2.input        = "infile:file, annofiles:files"
	pHeatmap2.output       = [
		"outfile:file:{{i.infile | fn2}}.heatmap/{{i.infile | fn2}}.heatmap.png",
		"outdir:dir:{{i.infile | fn2}}.heatmap"
	]
	pHeatmap2.args.devpars  = Box(res = 300, height = 2000, width = 2000)
	pHeatmap2.args.draw     = Box()
	pHeatmap2.args.params   = Box(heatmap_legend_param = Box())
	pHeatmap2.args.anopts   = Box(rnames = True, cnames = True)
	pHeatmap2.args.inopts   = Box(rnames = True, cnames = True)
	pHeatmap2.args.seed     = None
	pHeatmap2.args.saveinfo = True
	pHeatmap2.args.helper   = ''
	pHeatmap2.envs.rimport  = rimport
	pHeatmap2.lang          = params.Rscript.value
	pHeatmap2.script        = "file:scripts/plot/pHeatmap2.r"
	return pHeatmap2

@procfactory
def _pScatterCompare():
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
	return pScatterCompare

@procfactory
def _pROC():
	"""
	@input:
		infile: The input matrix file.
			- Col0: rownames if `args.inopts.rnames` is True
			- Col1: label (0, 1 class)
			- Col2: prediction values from model1
			- [Col3: prediction values from model2]
			- [...]
	@output:
		outdir: The output directory
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
	return Box(
		desc   = 'ROC curves and AUCs.',
		lang   = params.Rscript.value,
		input  = 'infile:file',
		output = [
			'outfile:file:{{i.infile | stem}}.roc/{{i.infile | stem}}.roc.png',
			'outdir:dir:{{i.infile | stem}}.roc'
		],
		args = Box(
			inopts  = Box(rnames = True, cnames = True),
			params  = Box(bestCut = True, showAUC = True),
			ggs     = Box({
				# show legend at bottom right corner
				'theme#auc': {
					'legend.position': [1, 0], 'legend.justification': [1, 0],
					'panel.border': 'r:element_rect(colour="black", fill=NA, size=.2)'}
			}),
			devpars = Box(res = 300, height = 2000, width = 2000)
		)
	)

@procfactory
def _pVenn():
	"""
	@description:
		Venn/UpsetR plots. Requires:
			- [`r-VennDiagram`](https://www.rdocumentation.org/packages/VennDiagram) and
			- [`r-UpSetR`](https://www.rdocumentation.org/packages/UpSetR)
	@input:
		`infile:file`: The input matrix, could be two formats:
			- `args.intype == "raw"`:
				```
				category1	category2	category3
				e1	e2	e2
				e2	e3	e4
				... ...
				```
			- `args.intype == "computed"`:
				```
					category1	category2	category3
				[e1]	0	1	1
				[e2]	0	0	1
				...
				[eN]	1	0	0
				```
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
		`intype`  : Type of input file. See `i.infile`. Default: `raw`
		`inopts`  : options to read the input file.
			- `rnames`  : Whether input file has rownames. Default: False
		`params`  : Other params for `venn.diagram` or `upset`. Default: {}
		`devpars` : The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`
	"""
	pVenn              = Proc(desc = 'Venn plots.')
	pVenn.input        = "infile:file, metafile:file"
	pVenn.output       = "outfile:file:{{i.infile | fn}}.venn.png"
	pVenn.args.tool    = 'auto' # upsetr or auto: < = 3 venn, else upsetr
	pVenn.args.inopts  = Box(rnames  = False, cnames = True)
	pVenn.args.intype  = 'raw' # computed
	pVenn.args.params  = Box()
	pVenn.args.devpars = Box(res = 300, height = 2000, width = 2000)
	pVenn.envs.rimport = rimport
	pVenn.lang         = params.Rscript.value
	pVenn.script       = "file:scripts/plot/pVenn.r"
	return pVenn

@procfactory
def _pVenn2():
	"""
	@name:
		pVenn2
	@description:
		Venn plots using individual input files, each of which contains the elements of the category.
	@input:
		`infiles:files`: The input files, each one is a category containing the elements.
			- If it has column name, then it will be used as category name, otherwise
			- the filename (without extension) will be used.
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
		`outfile:file`: The output file, Default: `{{i.infiles | fs2name}}.venn.png`
	@args:
		`tool`    : Which tools to use. Default: auto (venn, upsetr, auto(n<=3: venn, otherwise upsetr))
		`inopts`  : options to read the input file. Default: `Box(rnames = False, cnames = False)`
		`params`  : Other params for `venn.diagram` or `upset`. Default: `Box()`
		`devpars` : The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`
	"""
	pVenn2              = Proc(desc = 'Venn plots using individual input files')
	pVenn2.input        = "infiles:files, metafile:file"
	pVenn2.output       = "outfile:file:{{i.infiles | fs2name}}.venn.png"
	pVenn2.args.tool    = 'auto'
	pVenn2.args.inopts  = Box(rnames = False, cnames = False)
	pVenn2.args.params  = Box()
	pVenn2.args.devpars = Box(res = 300, height = 2000, width = 2000)
	pVenn2.envs.rimport = rimport
	pVenn2.envs.fs2name = fs2name
	pVenn2.lang         = params.Rscript.value
	pVenn2.script       = "file:scripts/plot/pVenn2.r"
	return pVenn2

@procfactory
def _pPie():
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
	return pPie

@procfactory
def _pManhattan():
	"""
	@name:
		pManhattan
	@description:
		Manhattan plot.
	@input:
		`infile:file`: The input file. First 6 columns should be BED6, and then column:
			- 7: The raw pvalue.
			- 8: The x-axis labels for the records. Optional. If not provided, chromosomes will be used.
			- This file has to be sorted by coordinates.
			- For example:
				```
				chr19	45604163	45604163	rs7255060	0	+	3.238E-03	+200K
				chr19	45595277	45595277	rs10417101	0	+	3.870E-03	+200K
				chr19	45394336	45394336	rs71352238	0	+	6.440E-03	-50K
				chr19	45615857	45615857	rs6509194	0	+	1.298E-02	+250K
				chr19	45594170	45594170	rs3178166	0	+	3.617E-02	+200K
				chr19	45361574	45361574	rs3852856	0	+	2.070E-02	-100K
				chr19	45220205	45220205	rs57090948	0	+	4.384E-02	-200K
				chr19	45396219	45396219	rs157582	0	+	9.472E-03	-50K
				chr19	45210634	45210634	rs10421830	0	+	1.375E-02	-250K
				chr19	45228502	45228502	rs10422350	0	+	4.121E-02	-200K
				```
		`hifile:file`: The file with the record names (one per line) to highlight in the plot.
	@output:
		`outfile:file`: The plot. Default: `{{i.infile | fn}}.manht.png`
	@args:
		`inopts` : Options to read the input file. Default: `Box(cnames = False, rnames = False)`
		`hilabel`: Show the labels of the highlight points. Default: `True`
		`ggs`    : Extra expressions for ggplot.
		`devpars`: The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`
		`gsize`  : The genome sizes file. Default: `None`
			- If `None`, will infer from the snp coordinates
	"""
	pManhattan              = Proc(desc = 'Manhattan plot.')
	pManhattan.input        = 'infile:file, hifile:file'
	pManhattan.output       = 'outfile:file:{{i.infile | fn}}.manht.png'
	pManhattan.args.inopts  = Box(cnames = False, rnames = False)
	pManhattan.args.hilabel = True
	pManhattan.args.devpars = Box(res = 300, height = 2000, width = 2000)
	pManhattan.args.ggs     = Box()
	pManhattan.args.gsize   = None # params.gsize.value
	pManhattan.envs.rimport = rimport
	pManhattan.lang         = params.Rscript.value
	pManhattan.script       = "file:scripts/plot/pManhattan.r"
	return pManhattan

@procfactory
def _pQQ():
	"""
	@name:
		pQQ
	"""
	pQQ              = Proc(desc = 'Q-Q plot')
	pQQ.input        = 'infile:file'
	pQQ.output       = 'outfile:file:{{i.infile | fn}}.qq.png'
	pQQ.args.inopts  = Box(cnames = True, rnames = False)
	pQQ.args.devpars = Box(res = 300, height = 2000, width = 2000)
	pQQ.args.ggs     = Box()
	pQQ.args.params  = Box()
	pQQ.envs.rimport = rimport
	pQQ.lang         = params.Rscript.value
	pQQ.script       = "file:scripts/plot/pQQ.r"
	return pQQ

@procfactory
def _pPairs():
	"""
	@name:
		pPairs
	"""
	pPairs              = Proc(desc = 'Plot pairs with ggpairs from GGally')
	pPairs.input        = 'infile:file'
	pPairs.output       = 'outfile:file:{{i.infile | fn}}.pairs.png'
	pPairs.args.inopts  = Box(cnames = True, rnames = True)
	pPairs.args.devpars = Box(res = 300, height = 400, width = 400) # for each cell
	pPairs.args.ggs     = Box() # geom_* not available but theme available
	pPairs.args.params  = Box(upper = Box(continuous = "density")) # other parameters for ggpairs
	pPairs.envs.rimport = rimport
	pPairs.lang         = params.Rscript.value
	pPairs.script       = "file:scripts/plot/pPairs.r"
	return pPairs

@procfactory
def _pVolcano():
	"""
	@name:
		pVolcano
	"""
	pVolcano               = Proc(desc = 'Do volcano plot.')
	pVolcano.input         = 'infile:file'
	pVolcano.output        = 'outfile:file:{{i.infile | fn}}.volcano.png'
	pVolcano.args.fccut    = 2
	pVolcano.args.pvalcut  = .05
	pVolcano.args.usepval  = False
	pVolcano.args.hilights = []
	pVolcano.args.devpars  = Box(res = 300, height = 2000, width = 2000)
	pVolcano.args.ggs      = Box()
	pVolcano.envs.rimport  = rimport
	pVolcano.lang          = params.Rscript.value
	pVolcano.script        = "file:scripts/plot/pVolcano.r"
	return pVolcano

