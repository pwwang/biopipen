{{rimport}}('__init__.r', 'plot.r')
library(ComplexHeatmap)

infile    = {{i.infile | R}}
outfile   = {{o.outfile | R}}
inopts    = {{args.inopts | R}}
devpars   = {{args.devpars | R}}
anopts    = {{args.anopts | R}}
annofiles = {{i.annofiles | R}}
seed      = {{args.seed | R}}
set.seed(seed)

data = read.table.inopts(infile, inopts)
annos = list()
for (annofile in annofiles) {
	annos = c(annos, list(read.table.inopts(annofile, anopts)))
}
if (length(annos) == 1) {
	annos = annos[[1]]
}
eval(parse(text = {{args.helper | repr}}))
params    = {{args.params | R}}
drawps    = {{args.draw | R}}

plot.heatmap2(data, outfile, params = params, draw = drawps, devpars = devpars)
