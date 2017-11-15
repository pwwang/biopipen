from pyppl import Proc, Box
from .utils import plot
from . import params

"""
Generate plot using given data
"""

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
	`header`:    `header` parameter for `read.table`, default: True
	`rownames`:  `row.names` parameter for `read.table`, default: 1
	`params`:    Other parameters for `boxplot`, default: ""
"""
pBoxplot                     = Proc(desc = 'Do boxplot.')
pBoxplot.input               = "datafile:file"
pBoxplot.output              = "outpng:file:{{in.datafile | fn}}.boxplot.png"
pBoxplot.lang                = "Rscript"
pBoxplot.args.header         = True
pBoxplot.args.rownames       = 1
pBoxplot.args.params         = Box()
pBoxplot.tplenvs.plotBoxplot = plot.boxplot.r
pBoxplot.script              = """
{{plotBoxplot}}
data = read.table ("{{in.datafile}}", sep="\\t", header={{args.header | Rbool}}, row.names={{args.rownames}}, check.names=F)
plotBoxplot(data, {{out.outpng | quote}}, {{args.params | Rlist}})
"""

"""
@name:
	pScatterPlot
@description:
	Scatter plots with more information
@input:
	`infile:file`: The input file.
	- Format:
	```
		X	Y	Size	Color
	A	1	1	1	1
	B	2	2	2	2
	```
	- Column 3,4 can be omitted
@output:
	`outfile:file`: The plot
@args:
	`colfunc`: The functions to generate colors. Default: `heat.colors`
		- Available: rainbow, heat.colors, terrain.colors, topo.colors, and cm.colors.
	`type`:    The type of the symbols. Default: `circles`
		- Available: circles, squares, rectangles, stars, thermometers, boxplots.
	`inches`:  Scale the largest symbol to this size. Default: 1/3
	`data`:    The columns for render the symbols. Default: 0 (a simple dot plot)
		- circles:       3 (radii)
		- squares:       3 (length of sides)
		- rectangles:    3:4 (widths and heights)
		- starts:        3:? (?>5, a matrix with three or more columns giving the lengths of the rays from the center of the stars.)
		- thermometers:  3:? (?=5|6, The first two columns give the width and height of the thermometer symbols. If there are three columns, the third is taken as a proportion: the thermometers are filled (using colour fg) from their base to this proportion of their height. If there are four columns, the third and fourth columns are taken as proportions and the thermometers are filled between these two proportions of their heights. The part of the box not filled in fg will be filled in the background colour (default transparent) given by bg.) 
		- boxplots:      3:7 (a matrix with five columns. The first two columns give the width and height of the boxes, the next two columns give the lengths of the lower and upper whiskers and the fifth the proportion (with a warning if not in [0,1]) of the way up the box that the median line is drawn.)
	`main`:   The title of the plot. Default: NULL (the file name)
	`xlab`:   The labels for x axis. Default: colnames(mat)[1]
	`ylab`:   The labels for y axis. Default: colnames(mat)[2]
	`text`:   Whether to show the text of the symbols (rownames). Default: TRUE
"""
pScatterPlot = Proc(desc = 'Scatter plots.')
pScatterPlot.input        = "infile:file"
pScatterPlot.output       = "outfile:file:{{infile | fn}}.scatter.png"
pScatterPlot.args.symbolsParams = {
	'inches':   0.33,
	'bg':       '#000000',
	'circles':  'r:data[,ncol(data)]',
	'main':     'NULL',
	'xlab':     'NULL',
	'ylab':     'NULL'
}
pScatterPlot.args.textParams = {
	'cex': 0.6
}
pScatterPlot.args.text    = "TRUE"
pScatterPlot.args._plotSymbols = plot.symbols.r
pScatterPlot.args._plotText    = plot.text.r
pScatterPlot.lang         = "Rscript"
pScatterPlot.script       = """
{{args._plotSymbols}}
png (file = "{{outfile}}", res=300, width=2000, height=2000)
data = read.table ("{{infile}}", sep="\\t", header = T, row.names = 1, check.names = F)
cnames        = colnames(data)
rnames        = rownames(data)
symbolsParams = {{args.symbolsParams | Rlist}}
textParams    = {{args.textParams | Rlist}}
if (is.null(symbolsParams$main)) {
	symbolsParams$main = "{{infile | fn}}"
}
if (is.null(symbolsParams$xlab)) {
	xlab = cnames[1]
}
if (is.null(symbolsParams$ylab)) {
	ylab = cnames[2]
}
plotSymbols(data[,1], data[,2], symbolsParams)
if ({{args.text | Rbool}}) {
	{{args._plotText}}
	textParams$labels = rnames
	plotText(data[,1], data[,2], textParams)
}
dev.off()
"""

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
pHeatmap                     = Proc(desc = 'Plot heatmaps.')
pHeatmap.input               = "infile:file"
pHeatmap.output              = "outfile:file:{{in.infile | fn}}.heatmap.png"
pHeatmap.args.ggs            = Box()
pHeatmap.args.devpars        = Box({'res': 300, 'height': 2000, 'width': 2000})
pHeatmap.args.dendro         = Box({'dendro': True})
pHeatmap.args.header         = True
pHeatmap.args.rownames       = 1
pHeatmap.args.rows           = 'all' # top:100, both:50, bottom:100, random-both:50, random:100
pHeatmap.args.cols           = 'all' # top:100, both:50, bottom:100, random-both:50, random:100
pHeatmap.tplenvs.plotHeatmap = plot.heatmap.r
pHeatmap.lang                = params.Rscript.value
pHeatmap.script              = "file:scripts/plot/pHeatmap.r"

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
@output:
	`outfile:file`: The plot
@args:
	`tool`: Which tools to use. Default: auto (venn, upsetr, auto(n<=3: venn, otherwise upsetr))
	`rownames`: Whether input file has rownames. Default: False
	`vennParams`: Other params for `venn.diagram`. Default: {}
	`upsetParams`: Other params for `upset`. Default: {}
@requires:
	[`r-VennDiagram`](https://www.rdocumentation.org/packages/VennDiagram)
	[`r-UpSetR`](https://www.rdocumentation.org/packages/UpSetR)
"""
pVenn                  = Proc(desc = 'Venn plots.')
pVenn.input            = "infile:file"
pVenn.output           = "outfile:file:{{infile | fn}}.venn.png"
pVenn.args.tool        = 'auto' # upsetr or auto: <=3 venn, else upsetr
pVenn.args.rownames    = 'NULL'
pVenn.args.vennParams  = {}
pVenn.args.upsetParams = {}
pVenn.args._plotVenn   = plot.venn.r
pVenn.args._plotUpset  = plot.upset.r
pVenn.lang             = 'Rscript'
pVenn.script           = """
library(UpSetR)
mat         = read.table ("{{infile}}", sep="\\t", header = TRUE, row.names = {{args.rownames | R}}, check.names = F)
nc          = ncol(mat)
cnames      = colnames(mat)
vennParams  = {{args.vennParams | Rlist}}
upsetParams = {{args.upsetParams | Rlist}}
if (nc<=3 && ({{args.tool | quote}} == 'auto' || {{args.tool | quote}} == 'venn')) {
	{{args._plotVenn}}
	plotVenn (mat, {{outfile | quote}}, vennParams)
} else {
	{{args._plotUpset}}
	plotUpset(mat, {{outfile | quote}}, upsetParams)
}
"""
	