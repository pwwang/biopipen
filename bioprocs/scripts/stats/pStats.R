{{rimport}}('__init__.r', 'plot.r')
pdf(NULL) # preventing Rplot.pdf in cwd
library(ggpubr)
# usage: https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html

options(stringsAsFactors=FALSE)
infile  = {{i.infile | R}}
outdir  = {{o.outdir | R}}
inopts  = {{args.inopts | R}}
types   = {{args.types | R}}
groups  = {{args.groups | R}}
ignore  = {{args.ignore | R}}
devpars = {{args.devpars | R}}
devpars2 = devpars
devpars2$width = 2*devpars$width

indata = read.table.inopts(infile, inopts)
ninst  = nrow(indata)


do_continuous = function (feat, groups, outdir) {
	#histp = plot.histo(data, 'return')
	statfile = file.path(outdir, paste0('feature-stat.', feat, '.txt'))
	testfile = file.path(outdir, paste0('feature-test.', feat, '.txt'))
	densfile = file.path(outdir, paste0('feature-plot.', feat, '.png'))
	statdata = data.frame()
	fdata = as.numeric(indata[, feat])
	statdata = rbind(statdata, list(
		Group  = '__ALL__',
		N      = ninst,
		Max    = max(fdata, na.rm = TRUE),
		Min    = min(fdata, na.rm = TRUE),
		Mean   = mean(fdata, na.rm = TRUE),
		Median = median(fdata, na.rm = TRUE),
		Sd     = sd(fdata, na.rm = TRUE)
	))
	testdata = data.frame()
	densdata = data.frame(values = fdata, Group = '__ALL__')
	colnames(densdata)[1] = feat
	for (grup in groups) {
		if (grup == feat) {next}
		lvls = levels(factor(indata[, grup]))
		for (lvl in lvls) {
			lvldata = as.numeric(indata[indata[, grup] == lvl, feat])
			densdf  = data.frame(values = lvldata, Group = paste0(grup, '[', lvl ,']'))
			colnames(densdf)[1] = feat
			densdata = rbind(densdata, densdf)
			statdata = rbind(statdata, data.frame(
				Group  = paste0(grup, ' [', lvl, ']'),
				N      = length(lvldata),
				Max    = max(lvldata, na.rm = TRUE),
				Min    = min(lvldata, na.rm = TRUE),
				Mean   = mean(lvldata, na.rm = TRUE),
				Median = median(lvldata, na.rm = TRUE),
				Sd     = sd(lvldata, na.rm = TRUE),
				check.names = FALSE))
		}
		# t, wilcox, anova, Kruskal-Wallis
		t = NULL
		tryCatch({
			t = t.test(fdata ~ indata[, grup])
		}, error = function(e) {
			t = list(p.value = '-')
		})
		w = NULL
		tryCatch({
			w = wilcox.test(fdata ~ indata[, grup])
		}, error = function(e) {
			w = list(p.value = '-')
		})
		a = summary(aov(fdata ~ indata[, grup]))[[1]]
		k = kruskal.test(fdata ~ indata[, grup])

		testdata = rbind(testdata, data.frame(
			Group = paste0(grup, ' [', paste(lvls, collapse = ', ') ,']'),
			t.test = paste0('p=', ifelse(is.numeric(t$p.value), round(t$p.value, 3), t$p.value)),
			wilcox.test = paste0('p=', ifelse(is.numeric(w$p.value), round(w$p.value, 3), w$p.value)),
			ANOVA = paste0('p=', round(a[1,5], 3)),
			Kruskal.Wallis.test = paste0('p=', round(k$p.value, 3)),
		check.names = FALSE))
	}
	# align the violin and boxplot
	dodge = ggplot2::position_dodge(width = 0.4)
	densplot = plot.density(densdata, 'return', devpars = devpars)
	vioplot = plot.violin(densdata, 'return', devpars = devpars, ggs = list(geom_boxplot = list(
		position = dodge,
		width = .1
	)), params = list(position = dodge))

	p = ggarrange(densplot, vioplot, ncol = 2)
	save.plot(p, densfile, devpars2)
	write.table(statdata, statfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
	if (nrow(testdata)) {
		write.table(testdata, testfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
	}
}

.count.levels = function(fdata, flevels) {
	counts = as.list(table(fdata))
	cnames = names(counts)
	ret = list()
	for (lvl in flevels) {
		ret[[lvl]] = ifelse(lvl %in% cnames, counts[[lvl]], 0)
	}
	ret
}

do_categorical = function(feat, groups, outdir) {
	# contingency tables
	statfile = file.path(outdir, paste0('feature-stat.', feat, '.txt'))
	testfile = file.path(outdir, paste0('feature-test.', feat, '.txt'))
	densfile = file.path(outdir, paste0('feature-plot.', feat, '.png'))
	fdata    = indata[, feat]
	flevels  = levels(fdata)
	statdata = data.frame(Group = '__ALL__')
	statdata = cbind(statdata, as.data.frame(.count.levels(fdata, flevels)))

	testdata = data.frame()
	for (grup in groups) {
		if (grup == feat) {next}
		lvls = levels(factor(indata[, grup]))
		for (lvl in lvls) {
			lvldata = indata[indata[, grup] == lvl, feat]

			counts = c(list(Group = paste0(grup, '[', lvl ,']')), .count.levels(lvldata, flevels))
			statdata = rbind(statdata, as.data.frame(counts))
		}
		# chi-squre, fisher exact
		csq = chisq.test(table(fdata, indata[, grup]))
		fex = fisher.test(table(fdata, indata[, grup]))

		testdata = rbind(testdata, data.frame(
			Group = paste0(grup, ' [', paste(lvls, collapse = ', ') ,']'),
			chi.squared.test = paste0('p=', round(csq$p.value, 3)),
			fisher.exact.test = paste0('p=', round(fex$p.value, 3)),
		check.names = FALSE))
	}

	plot.col(indata[, c(feat, grup), drop = FALSE], densfile, x = 2,
		params = list(ggplot2::aes_string(fill = feat)))
	write.table(statdata, statfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
	if (nrow(testdata)) {
		write.table(testdata, testfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
	}
}

features = make.names(colnames(indata))
colnames(indata) = features
features = setdiff(features, make.names(ignore))


# top 10 records of the data
write.table(
	indata[1:10, features, drop=FALSE],
	file.path(outdir, {{i.infile | stem | @append: ".top10.txt" | quote}}),
	row.names = inopts$rnames,
	col.names = inopts$cnames,
	sep = inopts$delimit,
	quote = FALSE)

for (grup in groups) {
	indata[, grup] = as.factor(as.character(indata[, grup]))
}
for (feat in features) {
	ftype = list.get(types, feat, 'auto')
	if (feat %in% groups) {
		ftype = 'categorical'
	} else if (ftype == 'auto') {
		lvls  = levels(factor(indata[, feat]))
		if (length(lvls) < 10 && length(lvls) < ninst * .8) {
			ftype = 'categorical'
		} else {
			lvls = as.numeric(indata[, feat])
			if (sum(is.na(lvls)) >= ninst * .9) {
				ftype = 'categorical'
			} else {
				ftype = 'continuous'
			}
		}
	}
	ftype = ifelse(startsWith(ftype, 'cont'), 'continuous', 'categorical')
	if (ftype == 'categorical') {
		indata[, feat] = as.factor(as.character(indata[, feat]))
	}
	do.call(paste0('do_', ftype), list(feat = feat, groups = make.names(groups), outdir = outdir))
}


