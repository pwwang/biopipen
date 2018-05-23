require('ggplot2')

.script.file = function() {
	for (i in 1:length(sys.frames())) {
		ofile = sys.frame(i)$ofile
		if (!is.null(ofile)) {
			return (ofile)
		}
	}
}
source(file.path(dirname(.script.file()), '__init__.r'))

# Update an aes in alist (named "mapping") for ggplot2 functions
# namesbefore = c('data')
update.aes = function(alist, newaes, namesbefore = c()) {
	alnames = names(alist)
	if ('mapping' %in% alnames) {
		alist$mapping = update.list(alist$mapping, newaes)
	} else {
		delidx = NULL
		if (length(namesbefore) > 0) {
			for (i in 1:length(namesbefore)) {
				if (namesbefore[i] %in% alnames)
					delidx = c(delidx, i)
			}
		}
		if (!is.null(delidx)) namesbefore = namesbefore[-delidx]
		pos = length(namesbefore) + 1
		if (is.null(alnames)) {
			if (length(alist) < pos) alist$mapping = newaes
			else alist[[pos]] = update.list(alist[[pos]], newaes)
		} else {
			x = 0
			if (length(alnames) > 0) {
				for (i in 1:length(alnames)) {
					if (alnames[i] != "") next
					x = x + 1
					if (x == pos) {
						alist[[i]] = update.list(alist[[i]], newaes)
						break
					}
				}
			}
			if (x == 0) {
				alist$mapping = newaes
			}
		}
	}
	return (alist)
}

apply.ggs = function(p, ggs) {
	if (is.null(ggs) || length(ggs) == 0)
		return(p)
	funcs = names(ggs)
	for (i in 1:length(funcs)) {
		p = p + do.call(funcs[i], ggs[[i]])
	}
	return (p)
}

plot.no = function(data, plotfile, params = list(), ggs = list(), devpars = list(res = 300, width = 2000, height = 2000)) {
	do.call(png, c(list(filename=plotfile), devpars))
	params$data = as.data.frame(data)
	p = do.call(ggplot, params)
	print(apply.ggs(p, ggs))
	dev.off()
}

plot.xy = function(data, plotfile, x = 1, y = 2, ggs = list(), devpars = list(res = 300, width = 2000, height = 2000)) {
	cnames = colnames(data)
	cnames = make.names(cnames)
	colnames(data) = cnames
	if (is.numeric(x)) x = cnames[x]
	if (is.numeric(y)) y = cnames[y]
	params = list(mapping = aes_string(x = x, y = y))
	plot.no(data, plotfile, params, ggs, devpars)
}

plot.x = function(data, plotfile, x = 1, ggs = list(), devpars = list(res = 300, width = 2000, height = 2000)) {
	cnames = colnames(data)
	cnames = make.names(cnames)
	colnames(data) = cnames
	if (is.numeric(x)) x = cnames[x]
	params = list(mapping = aes_string(x = x))
	plot.no(data, plotfile, params, ggs, devpars)
}

plot.stack = function(data, plotfile, x = 'ind', y = 'values', ggs = list(), devpars = list(res = 300, width = 2000, height = 2000)) {
	data = stack(as.data.frame(data))
	plot.xy(data, plotfile, x = x, y = y, ggs = ggs, devpars = devpars)
}

plot.roc = function(data, plotfile, params = list(returnAUC = T, showAUC = T, combine = T, labels = F), ggs = list(), devpars = list(res = 300, width = 2000, height = 2000)) {
	require('plotROC')

	cnames = colnames(data)
	if (is.null(cnames)) {
		cnames = c('D', paste0('M', 1:(ncol(data)-1)))
		colnames(data) = cnames
	}

	if (length(cnames) < 2) {
		stop('Only 1 column found in data!')
	} else if (length(cnames) == 2) {
		D = cnames[1]
		M = cnames[2]
		data$name = rep(M, nrow(data))
	} else {
		data = melt_roc(data, 1, 2:length(cnames))
		D = 'D'
		M = 'M'
	}

	returnAUC = params$returnAUC
	if (is.null(returnAUC)) returnAUC = T

	showAUC   = params$showAUC
	if (is.null(showAUC)) showAUC = T

	combine   = params$combine
	if (is.null(combine)) combine = T

	params$returnAUC = NULL
	params$showAUC   = NULL
	params$combine   = NULL

	params = update.aes(params, aes_string(d = D, m = M, color = 'name'))
	if (combine) {
		do.call(png, c(list(filename=plotfile), devpars))

		p = ggplot(data) + do.call(geom_roc, params)
		if (returnAUC || showAUC) {
			aucs = matrix(round(calc_auc(p)$AUC, 3), ncol = 1, nrow = length(cnames) - 1)
			rownames(aucs)  = cnames[-1]
			colnames(aucs)  = 'AUC'
		}
		if (showAUC) {
			if (length(cnames) > 2) {
				auclabels = NULL
				for (i in 1:length(cnames[-1])) {
					auclabels = c(auclabels, paste(cnames[-1][i], aucs[i, 1], sep = ': '))
				}
				p = p + scale_color_discrete(name = "AUC", breaks = cnames[-1], labels = auclabels)
			} else {
				p = p + annotate("text", x = .9, y = .1, label = paste('AUC', aucs, sep = ': '))
				p = p + scale_color_discrete(guide = F)
			}
		}
		print(apply.ggs(p, ggs))
		dev.off()
		if (returnAUC) return(aucs)
	} else {
		aucs = matrix(NA, ncol = 1, nrow = length(cnames) - 1)
		rownames(aucs)  = cnames[-1]
		colnames(aucs)  = 'AUC'
		for (cname in cnames[-1]) {
			prefix = tools::file_path_sans_ext(plotfile)
			do.call(png, c(list(filename = paste0(prefix, '-', cname, '.png')), devpars))
			p = ggplot(data[which(data$name == cname), , drop=F]) + do.call(geom_roc, params)
			if (returnAUC || showAUC) {
				aucs[cname, 1] = round(calc_auc(p)$AUC, 3)
			}
			if (showAUC) {
				p = p + annotate("text", x = .9, y = .1, label = paste('AUC', aucs[cname, 1], sep = ': '))
				p = p + scale_color_discrete(guide = F)
			}
			print(apply.ggs(p, ggs))
			dev.off()
		}
		if (returnAUC) return(aucs)
	}
}

plot.scatter = function(data, plotfile, x = 1, y = 2, params = list(), ggs = list(), devpars = list(res = 300, width = 2000, height = 2000)) {
	ggs = c(list(geom_point = params), ggs)
	plot.xy(data, plotfile, x, y, ggs, devpars)
}

# alias
plot.points = plot.scatter

plot.boxplot = function(data, plotfile, x = 1, y = 2, stack = F, params = list(), ggs = list(), devpars = list(res=300, width=2000, height=2000)) {
	if (!stack) {
		cnames = colnames(data)
		cnames = make.names(cnames)
		colnames(data) = cnames
		if (is.numeric(x)) {
			x = cnames[x]
		}
		params = update.aes(params, aes_string(group = x))

		ggs = c(
			list(geom_boxplot = params),
			list(theme = list(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))),
			ggs
		)
		plot.xy(data, plotfile, x, y, ggs, devpars)
	} else {
		ggs = c(
			list(geom_boxplot = params),
			list(theme = list(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))),
			ggs
		)
		plot.stack(data, plotfile, ggs = ggs, devpars = devpars)
	}
}

plot.heatmap = function(data, plotfile, params = list(dendro = T), ggs = list(), devpars = list(res=300, width=2000, height=2000)) {
	require('ggdendro')
	require('gtable')
	require('grid')

	do.call(png, c(list(filename=plotfile), devpars))

	theme_none = theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),
		axis.title.x = element_text(colour=NA),
		axis.title.y = element_blank(),
		axis.text.x = element_blank(),
		axis.text.y = element_blank(),
		axis.line = element_blank(),
		axis.ticks = element_blank()
	)

	ggplus = list(
		scale_fill_gradient2 = list(),
		theme = list(
			axis.ticks  = element_blank(),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(angle = 60, hjust = 1),
			legend.title = element_blank()
		)
	)

	dendro = params$dendro
	if (is.null(dendro)) dendro = T

	rnames = rownames(data)
	cnames = colnames(data)
	if (is.null(rnames)) {
		rnames = paste('ROW', 1:nrow(m), sep='')
		rownames(data) = rnames
	}
	if (is.null(cnames)) {
		cnames = paste('COL', 1:ncol(m), sep='')
		colnames(data) = cnames
	}
	# Dendrogram 1
	if (dendro == TRUE || dendro == 'both' || dendro == 'row') {
		ddrleft   = hclust(dist(data))
		plotdleft = ggplot(segment(dendro_data(as.dendrogram(ddrleft)))) +
		geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
		theme_none + theme(plot.margin = margin(r=-5, l = 10)) +
		coord_flip() +
		scale_y_reverse(expand = c(0, 0)) +
		scale_x_discrete(expand = c(0, 0.5))
		rowidx  = ddrleft$order
	} else {
		plotdleft = NULL
		rowidx    = 1:length(rnames)
		rowidx    = rev(rowidx)
	}

	# Dendrogram 2
	if (dendro == TRUE || dendro == 'both' || dendro == 'col') {
		ddrtop   = hclust(dist(t(data)))
		plotdtop = ggplot(segment(dendro_data(as.dendrogram(ddrtop)))) +
		geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
		theme_none + theme(plot.margin = margin(b=-15, t = 10)) +
		scale_x_discrete(expand = c(0, 0.5)) +
		scale_y_continuous(expand = c(0, 0))
		colidx  = ddrtop$order
	} else {
		plotdtop = NULL
		colidx   = 1:length(cnames)
	}

	mm        = stack(as.data.frame(data[rowidx, colidx]))
	mm$rnames = rnames[rowidx]
	plothm    = ggplot(mm, aes(ind, rnames)) + geom_tile(aes(fill=values)) +
	xlim(unique(as.vector(mm$ind))) + scale_y_discrete(position="right", limits=rnames)

	ggs = c(ggplus, ggs)
	plothm = apply.ggs(plothm, ggs)

	if (dendro == TRUE || dendro == 'both') {
		ghm = ggplotGrob(plothm)
		gd1 = ggplotGrob(plotdtop)
		gd2 = ggplotGrob(plotdleft)
		maxheights  = as.list(unit.pmax(ghm$heights, gd2$heights))
		ghm$heights = maxheights
		gd2$heights = maxheights
		maxwidths   = as.list(unit.pmax(ghm$widths, gd1$widths))
		ghm$widths  = maxwidths
		gd1$widths  = maxwidths

		g   = gtable(unit(c(.15,.85), "npc"), unit(c(.15,.85), "npc"))
		g   = gtable_add_grob(g, ghm, 2, 2)
		g   = gtable_add_grob(g, gd1, 1, 2)
		g   = gtable_add_grob(g, gd2, 2, 1)

		grid.newpage()
		grid.draw(g)
	} else if (dendro == 'col') {
		ghm = ggplotGrob(plothm)
		gd1 = ggplotGrob(plotdtop)

		maxwidths   = as.list(unit.pmax(ghm$widths, gd1$widths))
		ghm$widths  = maxwidths
		gd1$widths  = maxwidths

		g   = gtable(unit(1, "npc"), unit(c(.15,.85), "npc"))
		g   = gtable_add_grob(g, ghm, 2, 1)
		g   = gtable_add_grob(g, gd1, 1, 1)

		grid.newpage()
		grid.draw(g)
	} else if (dendro == 'row') {
		ghm = ggplotGrob(plothm)
		gd2 = ggplotGrob(plotdleft)
		maxheights  = as.list(unit.pmax(ghm$heights, gd2$heights))
		ghm$heights = maxheights
		gd2$heights = maxheights

		g   = gtable(unit(c(.15,.85), "npc"), unit(1, "npc"))
		g   = gtable_add_grob(g, ghm, 1, 2)
		g   = gtable_add_grob(g, gd2, 1, 1)

		grid.newpage()
		grid.draw(g)
	} else {
		print(plothm)
	}
	dev.off()
}

plot.histo = function(data, plotfile, x = 1, params = list(), ggs = list(), devpars = list(res=300, width=2000, height=2000)) {
	ggs = c(list(geom_histogram = params), ggs)
	plot.x(data, plotfile, x, ggs, devpars)
}

plot.freqpoly = function(data, plotfile, x = 1, params = list(), ggs = list(), devpars = list(res=300, width=2000, height=2000)) {
	ggs = c(list(geom_freqpoly = params), ggs)
	plot.x(data, plotfile, x, ggs, devpars)
}

plot.maplot = function(data, plotfile, ggs = list(), devpars = list(res=300, width=2000, height=2000)) {
	data = as.data.frame(data)
	cnames = colnames(data)
	A      = if ("A" %in% cnames) data$A else data[, 1]
	M      = if ("M" %in% cnames) data$M else data[, 2]
	thres  = if ("threshold" %in% cnames) data$threshold else if (length(cnames) > 2) data[, 3] else NULL
	data = data.frame(A, M, thres)

	ggs = c(list(
		theme       = list(legend.position = 'none'),
		geom_hline  = list(color = 'blue3', yintercept = 0, linetype = 'dashed'),
		stat_smooth = list(se = F, method = 'loess', color = 'red3')
	), ggs)

	params = list(size = 1.5, alpha = 1/5)
	if (!is.null(thres)) {
		params$mapping = aes(color = factor(thres))
	}
	plot.scatter(data, plotfile, x = 'A', y = 'M', params = params, ggs = ggs, devpars = devpars)
}

plot.venn = function(data, plotfile, params = list(), devpars = list(res=300, width=2000, height=2000)) {
	library(VennDiagram)
	rnames = rownames(data)
	if (is.null(rnames)) {
		rnames = paste0('ROW', 1:nrow(data))
		rownames(data) = rnames
	}
	x = list()
	for (cname in colnames(data)) {
		x[[cname]] = rownames(data[which(data[, cname] == 1), cname, drop=F])
	}
	default.params = list(x = x, filename = plotfile, height = devpars$height, width = devpars$width, resolution = devpars$res, imagetype = "png", alpha = .5, fill = rainbow(ncol(data)))
	params = update.list(default.params, params)
	do.call(venn.diagram, params)
}

plot.upset = function(data, plotfile, params = list(), devpars = list(res=300, width=2000, height=2000)) {
	library(UpSetR)
	do.call(png, c(list(filename = plotfile), devpars))
	default.params = list(data = data, sets = colnames(data))
	params = update.list(default.params, params)
	p = do.call(upset, params)
	print(p)
	dev.off()
}

plot.pie = function(data, plotfile, ggs = list(), devpars = list(res=300, width=2000, height=2000)) {
	percent    = function(x) paste0(format(round(x*100, 1), nsmall = 1), "%")
	N          = nrow(data)
	Group      = colnames(data)
	oValue     = if (N == 1) data[1,] else colSums(data)
	oValue     = as.vector(as.matrix(oValue))
	GroupName  = paste0(Group, " (", oValue, ")")
	Value      = 100*oValue/sum(oValue)
	Labels     = percent(Value / 100)
	data       = data.frame (Group = Group, Value = Value)
	data$Group = factor(data$Group, levels = data$Group)
	ytext      = cumsum(rev(Value)) - rev(Value) / 2

	default.ggs = list(
		geom_bar            = list(aes(fill = Group), width = 1, stat = 'identity'),
		geom_text           = list(y = ytext, label = rev(Labels)),
		coord_polar         = list('y'),
		scale_fill_discrete = list(labels = GroupName, breaks = Group)
	)
	ggs = c(default.ggs, ggs)
	plot.xy(data, plotfile, x = '""', y = 'Value', ggs, devpars)
}

plot.volplot = function(data, plotfile, ggs = list(), devpars = list(res=300, width=2000, height=2000)) {
	data   = as.data.frame(data)
	cnames = names(data)
	logfc  = if ("logFC" %in% cnames) data$logFC else data[, 1]
	fdr    = if ("FDR" %in% cnames) data$FDR else data[, 2]
	fdr    = -log10(fdr)

	# cutoffs
	logfccut    = if ("logFCCut" %in% cnames) data$logFCCut[1] else 2
	fdrcut      = if ("FDRCut" %in% cnames) data$FDRCut[1] else 0.05
	fdrcutlabel = round(fdrcut, 3)
	fdrcut      = -log10(fdrcut)

	threshold = as.factor(abs(logfc) > logfccut & fdr > fdrcut)

	xm = min(max(abs(logfc)), 10)
	ym = min(max(fdr), 15)
	if (xm <= logfccut) logfccut = 1

	data = data.frame(logfc, fdr, threshold)

	ggs = c(list(
		geom_point = list(alpha = .4, size = 1.75, aes(color = threshold)),
		geom_hline = list(yintercept = fdrcut, linetype = 'dashed', color = 'blue3'),
		geom_vline = list(xintercept = c(-logfccut, logfccut), linetype = 'dashed', color = 'red3'),
		xlim       = list(c(-xm, xm)),
		ylim       = list(c(0, ym)),
		geom_text  = list(x = xm, y = fdrcut, label = paste('p', '=', fdrcutlabel), vjust = -1, hjust = 1, color = 'blue3'),
		geom_text  = list(x = -logfccut, y = ym, label = paste(-logfccut, 'fold'),  vjust = 1, hjust = -0.1, color='red3'),
		geom_text  = list(x = +logfccut, y = ym, label = paste0('+', logfccut, ' ', 'fold'),  vjust = 1, hjust = -0.1, color="red3"),
		theme      = list(legend.position = "none"),
		xlab       = list('log2 Fold Change'),
		ylab       = list('-log10(Pvalue)')
	), ggs)

	plot.xy(data, plotfile, x = 'logfc', y = 'fdr', ggs = ggs, devpars = devpars)
}
