library(reshape2)
{{rimport}}('__init__.r', 'plot.r')

infile       = {{i.infile | quote}}
outdir       = {{o.outdir | quote}}
outfile      = {{o.outfile | quote}}
byrow        = {{args.byrow | R}}
inopts       = {{args.inopts | R}}
doplot       = {{args.plot | R}}
pval         = {{args.pval | R}}
method       = {{args.method | quote}}
outfmt       = {{args.outfmt | quote}}
params       = {{args.params | R}}
devpars      = {{args.devpars | R}}
plotfile     = {{o.outfile | prefix | @append: ".png" | quote}}
pvalfile     = {{o.outfile | prefix | @append: ".pval.txt" | quote}}
outpairfile  = {{o.outfile | prefix | @append: ".pairs.txt" | quote}}
outmatfile   = {{o.outfile | prefix | @append: ".mat.txt" | quote}}
pvalpairfile = {{o.outfile | prefix | @append: ".pval.pairs.txt" | quote}}
pvalmatfile  = {{o.outfile | prefix | @append: ".pval.mat.txt" | quote}}
groupfile    = {{i.groupfile | R}}
if (is.false(groupfile)) {
	groupfile = {{args.groupfile | R}}
}

data = read.table.inopts(infile, inopts)
if (byrow) data = t(data)

if (is.true(groupfile)) {
	if (file.exists(groupfile)) {
		groupdata = read.table(groupfile, row.names = NULL, col.names = FALSE)
		groups = levels(as.factor(groupdata[, 2]))
		if (length(groups) != 2) {
			stop("Only 2 groups allowed!")
		}
		group1 = groupdata[groupdata[,2] == groups[1], 1, drop = TRUE]
		group2 = groupdata[groupdata[,2] == groups[2], 1, drop = TRUE]
	} else { # group variable assigned directly
		groups = unlist(strsplit(groupfile, ";", fixed = TRUE))
		if (length(groups) != 2) {
			stop("Only 2 groups allowed!")
		}
		group1 = unlist(strsplit(groups[1], ",", fixed = TRUE))
		group2 = unlist(strsplit(groups[2], ",", fixed = TRUE))
	}
} else {
	group1 = colnames(data)
	group2 = group1
}
corrmat = matrix(NA, ncol = length(group2), nrow = length(group1))
pvalmat = matrix(NA, ncol = length(group2), nrow = length(group1))
rownames(corrmat) = group1
colnames(corrmat) = group2
rownames(pvalmat) = group1
colnames(pvalmat) = group2

if (pval) {
	for (v1 in group1) {
		for (v2 in group2) {
			corr = cor.test(data[, v1], data[, v2], method = method)
			corrmat[v1, v2] = corr$estimate
			corrmat[v2, v1] = corr$estimate
			pvalmat[v1, v2] = corr$p.value
			pvalmat[v2, v1] = corr$p.value
		}
	}
	write.table(pvalmat, pvalmatfile, row.names = T, col.names = T, quote = F, sep = "\t")
	ppairs = melt(pvalmat, na.rm = T)
	write.table(ppairs, pvalpairfile, row.names = F, col.names = F, quote = F, sep = "\t")
} else {
	corrmat = cor(data[, group1, drop = FALSE], data[, group2, drop = FALSE], method = method, use = "na.or.complete")
}
write.table(corrmat, outmatfile, row.names = T, col.names = T, quote = F, sep = "\t")
cpairs = melt(corrmat, na.rm = T)
write.table(cpairs, outpairfile, row.names = F, col.names = F, quote = F, sep = "\t")

if (outfmt == 'matrix') {
	file.symlink(outmatfile, outfile)
	if(pval) file.symlink(pvalmatfile, pvalfile)
} else {
	file.symlink(outpairfile, outfile)
	if(pval) file.symlink(pvalpairfile, pvalfile)
}

if (doplot) {
	col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
	plot.heatmap2(corrmat, plotfile = plotfile, params = c(list(
		col                  = col_fun,
		cluster_rows         = FALSE,
		cluster_columns      = FALSE,
		row_names_side       = "left",
		column_names_side    = "top",
		column_names_rot     = 60,
		heatmap_legend_param = list(title = "Correlation")), params
	), devpars = devpars)
}
