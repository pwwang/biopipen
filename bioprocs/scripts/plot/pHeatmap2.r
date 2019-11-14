{{rimport}}('__init__.r', 'plot.r')
library(ComplexHeatmap)

infile    = {{i.infile | R}}
outfile   = {{o.outfile | R}}
outdir    = {{o.outdir | R}}
hmfile    = {{o.outfile | prefix | @append: ".rds" | quote}}
rcfile    = {{o.outfile | prefix | @append: ".rowclusters.xls" | quote}}
ccfile    = {{o.outfile | prefix | @append: ".colclusters.xls" | quote}}
inopts    = {{args.inopts | R}}
devpars   = {{args.devpars | R}}
anopts    = {{args.anopts | R}}
saveinfo  = {{args.saveinfo | R}}
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

hm = plot.heatmap2(data, 'return', params = params, draw = drawps, devpars = devpars)
do.call(png, c(list(filename=outfile), devpars))
do.call(ComplexHeatmap::draw, c(list(hm), drawps))
dev.off()

if (saveinfo) {
	saveRDS(hm, hmfile)

	# export clusters
	ro = row_order(hm)
	# if no split, treat the whole as one cluster
	if (!is.list(ro)) ro = list(ro)
	if (is.null(names(ro))) {
		names(ro) = paste0('C', 1:length(ro))
	}
	rn_orig = rownames(data)
	if (is.null(rn_orig)) {
		rn_orig = 1:nrow(data)
	}
	rclines = c()
	for (clname in names(ro)) {
		rclines = c(
			rclines,
			paste("# Cluster:", clname, ', Size:', length(ro[[clname]])),
			paste(rn_orig[ ro[[clname]] ], collapse = ", ")
		)
	}
	rc_conn = file(rcfile)
	writeLines(rclines, rc_conn)
	close(rc_conn)

	co = column_order(hm)
	# if no split, treat the whole as one cluster
	if (!is.list(co)) co = list(co)
	if (is.null(names(co))) {
		names(co) = paste0('C', 1:length(co))
	}
	cn_orig = colnames(data)
	if (is.null(cn_orig)) {
		cn_orig = 1:nrow(data)
	}
	cclines = c()
	for (clname in names(co)) {
		cclines = c(
			cclines,
			paste("# Cluster:", clname, ', Size:', length(co[[clname]])),
			paste(cn_orig[ co[[clname]] ], collapse = ", ")
		)
	}
	cc_conn = file(ccfile)
	writeLines(cclines, cc_conn)
	close(cc_conn)
}