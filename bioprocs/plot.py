from pyppl import proc
from .utils import plot
"""
Generate plot using given data
"""

"""
@name:
	pBoxPlot
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
pBoxPlot = proc()
pBoxPlot.input  = "datafile:file"
pBoxPlot.output = "outpng:file:{{datafile | fn}}.png"
pBoxPlot.lang   = "Rscript"
pBoxPlot.args   = {"header": True, "rownames": 1, "params": ""}
pBoxPlot.script = """
png (file = "{{outpng}}", res=300, width=2000, height=2000)
data = read.table ("{{datafile}}", sep="\\t", header={{args.header | Rbool}}, row.names={{args.rownames}}, check.names=F)
boxplot (
	data{{args.params | lambda x: ", " + x if x else "" }}
)
dev.off()
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
pScatterPlot = proc (desc = 'Scatter plots.')
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

pHeatmap                   = proc (desc = 'Plot heatmaps.')
pHeatmap.input             = "infile:file"
pHeatmap.output            = "outfile:file:{{infile | fn}}.heatmap.png"
pHeatmap.args.params       = {
	'border_color': "#FFFFFF",
}
pHeatmap.args.transpose    = False
pHeatmap.args.header       = True
pHeatmap.args.rownames     = 1
pHeatmap.args.maxn         = 0
pHeatmap.args.subset       = 'random-all' # top, both, bottom, random-both
pHeatmap.args._plotHeatmap = plot.heatmap.r
pHeatmap.lang              = "Rscript"
pHeatmap.script            = """
{{args._plotHeatmap}}
params = {{args.params | Rlist}}
params$filename = {{outfile | quote}}
print ("Reading data ...")
mat = read.table ("{{infile}}", sep="\\t", header = {{args.header | Rbool}}, row.names = {{args.rownames}}, check.names = F)
if ({{args.transpose | Rbool}}) {
	print ("Transposing data ...")
	mat = t(mat)
}
maxn     = {{args.maxn}}
halfn    = as.integer(maxn/2)
nrows    = nrow(mat)
ncols    = ncol(mat)
nmax     = max(nrows, ncols)
if (nmax > maxn && maxn != 0) {
	if (nrows > ncols) {
		print ("Subsetting data by rows ...")
		if ("{{args.subset}}" == 'both') {
			mat = mat[order(rowSums(mat), decreasing = T),][c(1:halfn, (nrows-halfn):nrows), ]
		} else if ("{{args.subset}}" == 'top') {
			mat = mat[order(rowSums(mat), decreasing = T),][1:maxn, ]
		} else if ("{{args.subset}}" == 'bottom') {
			mat = mat[order(rowSums(mat), decreasing = T),][(nrows-maxn):nrows, ]
		} else if ("{{args.subset}}" == 'random-all') {
			mat = mat[sample(nrows, maxn), ]
		} else {
			gt0 = mat[rowSums(mat)>=0, ]
			lt0 = mat[rowSums(mat)<0, ]
			mat = rbind(gt0[sample(nrow(gt0), halfn), ], lt0[sample(nrow(lt0), halfn), ])
			rm(gt0, lt0)
		}
	} else {
		print ("Subsetting data by columns ...")
		if ("{{args.subset}}" == 'both') {
			mat = mat[, order(colSums(mat), decreasing = T)][, c(1:halfn, (ncols-halfn):ncols)]
		} else if ("{{args.subset}}" == 'top') {
			mat = mat[, order(colSums(mat), decreasing = T)][, 1:maxn]
		} else if ("{{args.subset}}" == 'bottom') {
			mat = mat[, order(colSums(mat), decreasing = T)][, (ncols-maxn):ncols]
		} else if ("{{args.subset}}" == 'random-all') {
			mat = mat[, sample(ncols, maxn)]
		} else {
			gt0 = mat[, colSums(mat)>=0]
			lt0 = mat[, colSums(mat)<0]
			mat = cbind(gt0[, sample(ncol(gt0), halfn)], lt0[, sample(ncol(lt0), halfn)])
			rm(gt0, lt0)
		}
	}
}
print ("Plotting data ...")
plotHeatmap(mat, params)
write.table(mat, "{{outfile | prefix}}.txt", quote=F)
rm(mat)
"""
	