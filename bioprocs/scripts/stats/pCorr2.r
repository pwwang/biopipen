library(methods)
library(data.table)
{{rimport}}('__init__.r', 'plot.r')

infile1 = {{i.infile1    | quote}}
infile2 = {{i.infile2    | quote}}
outdir  = {{o.outdir     | quote}}
outfile = {{o.outfile    | quote}}
dopval  = {{args.pval    | R}}
dofdr   = {{args.fdr     | R}}
inopts1 = {{args.inopts1 | R}}
inopts2 = {{args.inopts2 | R}}
doplot  = {{args.plot    | R}}
method  = {{args.method  | quote}}
outfmt  = {{args.outfmt  | quote}}
params  = {{args.params  | R}}
ggs     = {{args.ggs     | R}}
devpars = {{args.devpars | R}}

indata1 = read.table.inopts(infile1, inopts1)
indata2 = read.table.inopts(infile2, inopts2)
len1    = ncol(indata1)
len2    = ncol(indata2)
indata = merge(indata1, indata2)

rnames1 = rownames(indata1)
rnames2 = rownames(indata2)
if (is.null(rnames1))
	rnames1 = paste0('R1_', 1:nrow(indata1))
if (is.null(rnames2))
	rnames2 = paste0('R2_', 1:nrow(indata2))
rm(indata1, indata2)
saveResults = function (corr, rnames1, rnames2, outfmt, outfile) {
	if (outfmt == 'pairs') {
		rownames(corr) = apply(merge(rnames1, rnames2), 1, function(row) paste(row, collapse = "\t"))
		write.table(corr, outfile, col.names = F, row.names = T, quote = F, sep = "\t")
	} else {
		corr = cbind(merge(rnames1, rnames2), corr = corr)
		corr = dcast(corr, x ~ y, value.var = 'corr')
		rownames(corr) = corr[, 1]
		corr = corr[, -1, drop = F]
		write.table(corr, outfile, col.names = T, row.names = T, quote = F, sep = "\t")
	}
}

if (!dopval && dofdr == F) {
	corr = apply(indata, 1, function(row) cor(
		row[1:len1],
		row[len1+1:len2],
		method = method
	))
	corr = data.frame(corr = corr)
	saveResults(corr, rnames1, rnames2, outfmt, outfile)
} else {
	corr0 = apply(indata, 1, function(row) cor.test(
		row[1:len1],
		row[len1+1:len2],
		method = method
	))
	corr = data.frame(corr = sapply(corr0, function(x) x$estimate))
	saveResults(corr, rnames1, rnames2, outfmt, outfile)
	pvals = sapply(corr0, function(x) x$p.value)
	corrp = data.frame(pval = sapply(corr0, function(x) x$p.value))
	saveResults(corrp, rnames1, rnames2, outfmt, paste0(tools::file_path_sans_ext(outfile), '_pval.txt'))

	if (dofdr != F) {
		if (dofdr == T) dofdr = 'BH'
		qvals = p.adjust(pvals, method = dofdr)
		corrq = data.frame(qval = qvals)
		saveResults(corrp, rnames1, rnames2, outfmt, paste0(tools::file_path_sans_ext(outfile), '_qval.txt'))
	}
}

if (doplot) {
	outplot = paste0(tools::file_path_sans_ext(outfile), '.png')
	if (is.null(params$dendro))
		params$dendro = F
	corr = cbind(merge(rnames1, rnames2), corr = corr)
	corr = dcast(corr, x ~ y, value.var = 'corr')
	rownames(corr) = corr[, 1]
	corr = corr[, -1, drop = F]
	#corr[lower.tri(corr)] = 0
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