library(reshape2)
{{rimport}}('__init__.r', 'plot.r')

infile   = {{i.infile    | quote}}
outdir   = {{o.outdir   | quote}}
outfile  = {{o.outfile  | quote}}
byrow    = {{args.byrow   | R}}
inopts   = {{args.inopts  | R}}
doplot   = {{args.plot    | R}}
method   = {{args.method  | quote}}
outfmt   = {{args.outfmt  | quote}}
params   = {{args.params  | R}}
ggs      = {{args.ggs     | R}}
devpars  = {{args.devpars | R}}

inparams = list(
	file        = infile,
	header      = ifelse(is.null(inopts$cnames) || inopts$cnames, T, F),
	row.names   = ifelse(is.null(inopts$rnames) || inopts$rnames, 1, NULL),
	sep         = ifelse(is.null(inopts$delimit), "\t", inopts$delimit),
	check.names = F
)

data = do.call(read.table.nodup, inparams)

if (byrow) data = t(data)

corr = cor(data, method = method)
corr = round(corr, 3)

if (outfmt == 'matrix') {
	write.table(corr, outfile, row.names = T, col.names = T, quote = F, sep = "\t")
} else {
	cmat = melt(corr, na.rm = T)
	write.table(cmat, outfile, row.names = F, col.names = F, quote = F, sep = "\t")
}

if (doplot) {
	outplot = file.path(outdir, paste(basename(outdir), method, 'png', sep = '.'))
	if (is.null(params$dendro))
		params$dendro = F
	corr[lower.tri(corr)] = 0
	ggs = c(list(
		theme = list(
			legend.position      = c(0.03, 0.97),
			legend.justification = c("left", "top"),
			legend.title         = element_text(margin = margin(b = 5)),
			legend.box.just      = "right"
		),
		guides = list(
			fill = guide_legend(title = paste0(toupper(substring(method, 1, 1)), substring(method, 2)), reverse = T)
		)
	), ggs)
	plot.heatmap(corr, outplot, params, ggs, devpars)
}
