from pyppl import proc

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
	