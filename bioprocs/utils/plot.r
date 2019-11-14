require('ggplot2')
require('ggrepel')
pdf(NULL) # preventing Rplots.pdf

.script.file = function() {
	for (i in 1:length(sys.frames())) {
		ofile = sys.frame(i)$ofile
		if (!is.null(ofile)) {
			return (ofile)
		}
	}
}
source(file.path(dirname(.script.file()), '__init__.r'))

DEVPARS = list(height = 2000, width = 2000, res = 300)

.stack.data = function(data, colnames = c('values', 'ind')) {
	data = stack(as.data.frame(data))
	colnames(data) = colnames
	return(data)
}

.get.aes.string = function(data, ...) {
	# colnames(data) = c('A', 'B', 'C')
	# .get.aes.string(data, x=2, y=1)
	# => list(x='B', y='A')
	cols   = list(...)
	cnames = colnames(data)
	ret    = list()
	anames = names(cols)
	for (i in 1:length(anames)) {
		aname = anames[i]
		colidx = cols[[aname]]
		ret[[aname]] = ifelse(is.numeric(colidx), cnames[colidx], cols[[aname]])
	}
	return (ret)
}

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
		if (!is.null(ggs[[i]])) {
			p = p + do.call(funcs[i], ggs[[i]])
		}
	}
	return (p)
}

save.plot = function(p, plotfile, devpars = DEVPARS) {
	if (is.null(plotfile)) {
		print(p)
	} else if (plotfile == 'return') {
		return(p)
	} else if (endsWith(plotfile, '.png')) {
		do.call(png, c(list(filename = plotfile), devpars))
		print(p)
		dev.off()
	} else {
		stop('Only .png file supported to save plots for now.')
	}
}

plot.no = function(data, plotfile = NULL, params = list(), ggs = list(), devpars = DEVPARS) {
	params$data = as.data.frame(data)
	p = do.call(ggplot, params)
	p = apply.ggs(p, ggs)
	save.plot(p, plotfile, devpars)
}

plot.xy = function(data, plotfile = NULL, x = 1, y = 2, ggs = list(), devpars = DEVPARS) {
	cnames = colnames(data)
	cnames = make.names(cnames)
	colnames(data) = cnames
	if (is.numeric(x)) x = sprintf("`%s`", cnames[x])
	if (is.numeric(y)) y = sprintf("`%s`", cnames[y])
	params = list(mapping = aes_string(x = x, y = y))
	plot.no(data, plotfile, params, ggs, devpars)
}

plot.x = function(data, plotfile, x = 1, ggs = list(), devpars = DEVPARS) {
	cnames = colnames(data)
	cnames = make.names(cnames)
	colnames(data) = cnames
	if (is.numeric(x)) x = cnames[x]
	params = list(mapping = aes_string(x = x))
	plot.no(data, plotfile, params, ggs, devpars)
}

plot.stack = function(data, plotfile, x = 'ind', y = 'values', ggs = list(), devpars = DEVPARS) {
	data = stack(as.data.frame(data))
	plot.xy(data, plotfile, x = x, y = y, ggs = ggs, devpars = devpars)
}

plot.roc = function(
	data,
	plotfile = NULL,
	stacked = FALSE,
	params = list(returnTable = FALSE, showAUC = TRUE, bestCut = FALSE),
	ggs = list(), devpars = DEVPARS) {
	# plot the ROC curve
	# see: https://cran.r-project.org/web/packages/plotROC/vignettes/examples.html
	# @params:
	#	`data`: The data for plotting.
	#		- If stacked, then there should be only 3 columns: D, M, name
	#		- Else the 1st column should be D, and rest are `<M1>, <M2>, ...`
	#		- D must be binary.
	#	`plotfile`: The file to save the plot.
	#	`stacked` : Whether the data is stacked(melt). See `data`
	#	`params`  : The parameters for plotting.
	#		- `showAUC`  : Show AUC on the plot
	#		- `returnTable`: Whether return a table of details. If TRUE, to get plot, use ret$plot and table use ret$table
	require('pROC')

	if (stacked) {
		if (ncol(data) != 3)
			stop('Expect 3 columns (D, M, Group) in stacked data to plot ROC.')
		colnames(data) = c('D', 'M', 'Group')
	} else {
		ncols = ncol(data)
		if (ncols < 2) {
			stop('Expect at least 2 columns (D, M1, ..., Mn)')
		}
		data = cbind(rep(data[,1], ncols-1), stack(data[, 2:ncols, drop = FALSE]))
		colnames(data) = c('D', 'M', 'Group')
	}
	groups = levels(factor(data$Group))
	returnTable = params$returnTable
	if (is.null(returnTable)) returnTable = FALSE

	showAUC = params$showAUC
	if (is.null(showAUC)) showAUC = TRUE

	bestCut = params$bestCut
	if (is.null(bestCut)) bestCut = FALSE

	params$returnTable = NULL
	params$showAUC     = NULL
	params$bestCut     = NULL

	# get roc objects
	rocs = list()
	for (grup in groups) {
		rocs[[grup]] = roc(D ~ M, data = data, subset=(Group == grup))
	}

	# get table of a roc object
	roc.table = function(rocobj) {
		aucci = ci(rocobj)
		ret = data.frame(
			threshold   = rocobj$thresholds,
			sensitivity = rocobj$sensitivities,
			specificity = rocobj$specificities,
			auc         = as.vector(rocobj$auc),
			auc_95ci1   = aucci[1],
			auc_95ci2   = aucci[3])
		# checkout the best cut points
		bests = as.matrix(coords(rocobj, 'best'))
		cbind(ret, is.best = ret$threshold %in% bests[1,])
	}
	roc.tables = function(rocs) {
		ret = NULL
		for (name in names(rocs)) {
			rtable = roc.table(rocs[[name]])
			tmp = cbind(Group = name, rtable)
			ret = ifelse(is.null(ret), tmp, rbind(ret, tmp))
		}
		ret
	}

	ret       = list()
	ret$plot  = do.call(ggroc, c(list(rocs), params))
	ret$table = roc.tables(rocs)
	# add best cut
	if (bestCut) {
		bestdata = ret$table[ret$table$is.best,,drop=FALSE]
		ret$plot = ret$plot + geom_point(
			aes(x = specificity, y = sensitivity, color = Group),
			data = bestdata,
			inherit.aes = FALSE
		) + geom_segment(
			aes(x = specificity, xend = specificity, y = sensitivity, yend = -Inf, color = Group),
			data = bestdata,
			linetype = 'dashed',
			inherit.aes = FALSE
		) + geom_segment(
			aes(x = Inf, xend = specificity, y = sensitivity, yend = sensitivity, color = Group),
			data = bestdata,
			linetype = 'dashed',
			inherit.aes = FALSE
		) + geom_text_repel(
			aes(x = specificity, y = sensitivity,
				label = paste0('best (', round(specificity, 3), ', ', round(sensitivity, 3) ,')'),
				color = Group),
			data = bestdata
		)
	}
	# add AUC
	if (showAUC) {
		aucs = unique(ret$table[, c('Group', 'auc')])
		if (nrow(aucs) == 1) {
			ret$plot = ret$plot + annotate("text", x = .1, y = .05,
				label = paste('AUC', round(aucs[1,2], 3), sep = ' = '))
			ret$plot = ret$plot + scale_color_discrete(guide = F)
		} else {
			auclabels = apply(aucs, 1, function(row) {
				sprintf('AUC(%s) = %.3f', row[1], as.numeric(row[2]))
			})
			ret$plot = ret$plot + scale_color_discrete(name = "",
				breaks = groups,
				labels = as.vector(auclabels))
		}
	}

	ret$plot = apply.ggs(ret$plot, ggs)
	if (returnTable) {
		save.plot(ret$plot, plotfile, devpars)
		return(ret)
	}
	save.plot(ret$plot, plotfile, devpars)
}

plot.scatter = function(data, plotfile = NULL, x = 1, y = 2, params = list(), ggs = list(), devpars = DEVPARS) {
	ggs = c(list(geom_point = params), ggs)
	plot.xy(data, plotfile, x, y, ggs, devpars)
}
# alias
plot.points = plot.scatter

plot.col = function(
	data,
	plotfile = NULL,
	x = 2,
	y = 0,
	params = list(),
	ggs = list(), devpars = DEVPARS) {
	# y = 1 # y as counts
	# x = 2
	# ----------------------
	# Counts	Group	Fill
	# 10	Ga	Yes
	# 10	Gb	No
	#
	#    |
	# 10 | |-|  |-|     Fill
	#    | |x|  |.|     x: Yes
	#    |_|x|__|.|___  .: No
	#      Ga    Gb
	#
	# y.is.counts = FALSE
	# x = 1
	# y = 0 # counting Sub-groups
	# ---------------------------
	# Group	Fill
	# Ga	Yes
	# Ga	Yes
	# Ga	Yes
	# Gb	No
	# Gb	No
	# Gb	No
	#
	#    |
	# 3  | |-|  |-|     Fill
	#    | |x|  |.|     x: Yes
	#    |_|x|__|.|___  .: No
	#      Ga    Gb
	#
	if (y == 0) {
		data = cbind(data, Count = 1)
		y = ncol(data)
	}

	aes.string = .get.aes.string(data, x = x, y = y)
	ggs = c(
		list(geom_col = params),
		list(theme = list(axis.text.x = element_text(angle = 60, hjust = 1))),
		list(xlab = list(aes.string$x)),
		ggs
	)
	data[, x] = as.character(data[, x])
	plot.xy(data, plotfile, aes.string$x, aes.string$y, ggs, devpars)
}
# alias
plot.bar = plot.col

plot.boxplot = function(data, plotfile = NULL, x = 2, y = 1, stacked = TRUE, params = list(), ggs = list(), devpars = DEVPARS) {
	if (stacked) {
		cnames = colnames(data)
		cnames = make.names(cnames)
		colnames(data) = cnames
		if (is.numeric(x)) {
			x = cnames[x]
		}
		if (is.numeric(y)) {
			y = cnames[y]
		}
		#params = update.aes(params, aes_string(group = x))

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

plot.violin = function(data, plotfile = NULL, x = 2, y = 1, stacked = TRUE, params = list(), ggs = list(), devpars = DEVPARS) {
	if (!stacked) {
		data = .stack.data(data)
		x = 2
		y = 1
	}
	aes.string = .get.aes.string(data, x = x, y = y)
	ggs = c(
		list(geom_violin = params),
		list(theme = list(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))),
		ggs
	)
	plot.xy(data, plotfile, aes.string$x, aes.string$y, ggs, devpars)
}

plot.heatmap2 = function(
	data, plotfile = NULL, params = list(), draw = list(),
	devpars = DEVPARS) {
	library(ComplexHeatmap)

	params$matrix = data
	hm = do.call(Heatmap, params)

	if (is.null(plotfile)) {
		do.call(ComplexHeatmap::draw, c(list(hm), draw))
	} else if (plotfile == 'return') {
		return(hm)
	}  else {
		do.call(png, c(list(filename=plotfile), devpars))
		do.call(ComplexHeatmap::draw, c(list(hm), draw))
		dev.off()
	}
}

plot.heatmap = function(data, plotfile, params = list(dendro = T), ggs = list(), devpars = DEVPARS) {
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

plot.histo = function(data, plotfile = NULL, x = 1, params = list(), ggs = list(), devpars = DEVPARS) {
	ggs = c(list(geom_histogram = params), ggs)
	plot.x(data, plotfile, x, ggs, devpars)
}

plot.density = function(
	data,            # the data, either stacked or not
	# NULL: print the plot,
	# 'return': return the ggplot object,
	# 'filepath': save the plot to file
	plotfile = NULL,
	x = 1, # The data column for values, only for stacked data
	y = 2, # The data column for groups, only for stacked data
	stacked = TRUE,  # whether the data is stacked, if not will stack it
	params = list(), # The params for geom_density
	ggs = list(), devpars = DEVPARS) {

	if (!stacked) {
		data = .stack.data(data)
		x = 1
		y = 2
	}
	aes.string = .get.aes.string(data, x = x, y = y)
	params = c(params, list(mapping = aes_string(fill = aes.string$y), alpha = .3))
	ggs = c(list(geom_density = params), ggs)
	plot.x(data, plotfile, aes.string$x, ggs, devpars)
}


plot.freqpoly = function(data, plotfile, x = 1, params = list(), ggs = list(), devpars = DEVPARS) {
	ggs = c(list(geom_freqpoly = params), ggs)
	plot.x(data, plotfile, x, ggs, devpars)
}

plot.maplot = function(data, plotfile, threshold, ggs = list(), devpars = DEVPARS) {
	data = as.data.frame(data)
	cnames = colnames(data)
	A      = if ("A" %in% cnames) data$A else data[, 1]
	M      = if ("M" %in% cnames) data$M else data[, 2]
	thres  = threshold
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

plot.qq = function(data, plotfile = NULL, x = NULL, y = 1, params = list(), ggs = list(), devpars = DEVPARS) {
	data = as.data.frame(data)
	n    = nrow(data)
	q    = (1:n)/n
	if (is.null(x)) {
		data = cbind(data, Theoretical = qnorm(q))
		x = "Theoretical"
	}
	data[, x] = quantile(data[, x], q)
	data[, y] = quantile(data[, y], q)
	if (is.numeric(x)) x = sprintf("`%s`", colnames(data)[x])
	if (is.numeric(y)) y = sprintf("`%s`", colnames(data)[y])

	plot.scatter(data, plotfile, x = x, y = y, params = params, ggs = ggs, devpars = devpars)
}

# data:
# 		cat1	cat2	cat3
# item1	1		0		1
# item2	1		1		0
# item3	0		1		1
# item4	1		0		0
# item5	1		0		1
plot.venn = function(data, plotfile = NULL, params = list(), devpars = DEVPARS) {

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

plot.upset = function(data, plotfile = NULL, params = list(), devpars = DEVPARS) {
	library(UpSetR)
	default.params = list(data = data, nsets = ncol(data))
	params = update.list(default.params, params)
	p = do.call(upset, params)
	save.plot(p, plotfile, devpars)
}

plot.pie = function(data, plotfile, ggs = list(), devpars = DEVPARS) {
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

# Volcano plot
# data must be a data frame of log fold change and p/q value
# params is a list of:
#	- logfccut: The absolute log fold change cutoff
#	- pcut: The p/q value cutoff
#	- hilight: genes to highlight
#		- either the names (rownames) of the genes, or
#		- A number or a vector of two numbers, of genes to highlight for both up-/down genes, or
#		- a list of cutoffs for logfc and p/q value
#			- list(logfc = 2, p = .05), or if you want have different logfc cutoff for up-/down genes:
#			- list(logfc = c(-2, +2), p = .05), or
#
plot.volplot = function(data,
	plotfile,
	params  = list(logfccut = 2, pcut = 0.05, hilight = 5),
	ggs     = list(),
	devpars = DEVPARS) {

	data = as.data.frame(data)
	# let's see if it is p or q value
	is.pval = grepl('p', colnames(data)[2], fixed = TRUE)
	colnames(data) = ifelse(is.pval, c('logFC', 'log10.P.'), c('logFC', 'log10.FDR.'))
	data[, 2] = -log10(data[, 2])

	pcutlabel = round(params$pcut, 3)
	params$pcut = -log10(params$pcut)
	data$Group = apply(data, 1, function(row) {
		if (row[1]>=params$logfccut && row[2]>params$pcut) {
			'UP_SIG'
		} else if (row[1]<=-params$logfccut && row[2]>params$pcut) {
			'DOWN_SIG'
		} else {
			'INSIG'
		}
	})
	data$Group = factor(data$Group, levels = c('DOWN_SIG', 'UP_SIG', 'INSIG'))
	# parse Labels
	rnames = rownames(data)
	labeldata = NULL
	if (is.null(params$hilight)) {
		# pass
	} else if (!is.list(params$hilight) && !is.numeric(params$hilight)) {
		# gene names
		labeldata = data[params$hilight, 1:2]
		labeldata = cbind(labeldata, Label = params$hilight)
	} else if (!is.list(params$hilight) && is.numeric(params$hilight)) {
		# tops
		if (length(params$hilight) == 1) {
			params$hilight = c(params$hilight, params$hilight)
		}
		downtop   = params$hilight[1]
		uptop     = params$hilight[2]
		labeldata = data[order(data$logFC),,drop = FALSE]
		labeldata = labeldata[labeldata$Group != 'INSIG',,drop = FALSE]
		ngenes    = nrow(labeldata)
		labeldata = labeldata[c(1:downtop, (ngenes-uptop+1):ngenes),,drop = FALSE]
		labeldata$Label = rownames(labeldata)
	} else if (is.list(params$hilight)) {
		logfc = params$hilight$logfc
		if (is.null(logfc)) {
			logfc = c(-4, 4)
		} else if (length(logfc) == 1) {
			logfc = c(-logfc, logfc)
		}
		p         = params$hilight$p
		p         = -log10(ifelse(is.null(p), 0.01, p))
		labeldata = data[(data$logFC <= logfc[1] | data$logFC >= logfc[2]) & data[,2] > p,,drop = FALSE]
		labeldata$Label = rownames(labeldata)
	}

	vline = list(xintercept = c(-params$logfccut, params$logfccut),
		linetype = 'dashed', color = c('red3', 'blue3'))
	vtext1 = list(x = -params$logfccut, y = 0,
		label = paste(-params$logfccut, 'fold'), vjust = 1, hjust = 1.1, color='red3')
	vtext2 = list(x = +params$logfccut, y = 0,
		label = paste0('+', params$logfccut, ' ', 'fold'), vjust = 1, hjust = -0.1, color="blue3")
	hline = list(yintercept = params$pcut, linetype = 'dashed', color = 'black')
	htext = list(x = max(data$logFC), y = params$pcut,
		label = paste(ifelse(is.pval, 'P', 'FDR'), '=', pcutlabel),
		vjust = -1, hjust = 1, color = 'black')
	labels = NULL
	if (!is.null(labeldata)) {
		labels = list(aes_string(
			x = 'logFC', y = ifelse(is.pval, 'log10.P.', 'log10.FDR.'), label = 'Label', color = 'Group'
		), alpha = 1, data = labeldata, min.segment.length = unit(0, 'lines'), inherit.aes = FALSE)
	}
	xlim = max(max(data$logFC), abs(min(data$logFC)))
	ggs = c(list(
		geom_point         = list(alpha = .4, size = 1.75, aes(color = Group)),
		scale_color_manual = list(values=c("#F8766D", "#619CFF", "#999999")),
		geom_hline         = hline,
		geom_text          = htext,
		geom_vline         = vline,
		geom_text          = vtext1,
		geom_text          = vtext2,
		xlim               = list(c(-xlim, xlim)),
		theme_bw           = list(),
		theme              = list(legend.position ="none",
								panel.grid.major = element_blank(),
								panel.grid.minor = element_blank()),
		geom_label_repel = labels,
		xlab             = list('log2 Fold Change'),
		ylab             = list(ifelse(is.pval, '-log10(P)', '-log10(FDR)'))
	), ggs)

	plot.xy(data, plotfile,
		x = 'logFC', y = ifelse(is.pval, 'log10.P.', 'log10.FDR.'),
		ggs = ggs, devpars = devpars)
}

plot.man = function(data, plotfile = NULL, hilights = list(), hilabel = TRUE,
	gsize = NULL, ggs = list(), devpars = DEVPARS) {
	# manhattan plot
	# data is a data frame of
	# Chr, Pos, P[, Region]
	# Rownames should be snp name
	# The region is used to divide the plot, in case
	# all snps are on the same chromosome

	library(data.table)
	data = as.data.table(data)
	if (ncol(data) == 4) {
		colnames(data) = c('Snp', 'Chr', 'Pos', 'P')
		data$Region = data$Chr
	} else if (ncol(data) == 5) {
		colnames(data) = c('Snp', 'Chr', 'Pos', 'P', 'Region')
	} else {
		stop('Expect a data frame of 4 or 5 columns to do manhattan plot.')
	}
	# calculate the x axis for each snp, position on each chromsome should be accumulated.
	# min and max pos for each chr

	if (is.null(gsize)) {
		# Chr	minPos	maxPos	cumPos
		# chr1	910473	249162263	0
		# chr2	449487	242926381	249162263
		# chr3	281636	197837967
		chrlen = data[, .(minPos = min(Pos), maxPos = max(Pos)), by = Chr][, .(Chr, minPos, maxPos, cumPos = shift(cumsum(as.numeric(maxPos)), 1, fill = 0))]
		# add it back to data
	} else {
		chrs   = unique(unlist(data[, "Chr"]))
		gsize  = gsize[chrs,,drop = FALSE]
		gsize  = as.data.table(cbind(Chr = rownames(gsize), minPos = 0, maxPos = gsize[,1]))
		gsize  = gsize[Chr %in% chrs]
		chrlen = gsize[, .(Chr, minPos, maxPos, cumPos = shift(cumsum(as.numeric(maxPos)), 1, fill = 0))]
	}
	data = merge(data, chrlen, by = "Chr", sort = FALSE)

	data[, X:= Pos + cumPos][, Y:=-log10(P)]
	# get region centers
	rdata = data[, .(Region, centerPos = (min(X) + max(X))/2), by = Region]
	# fix the order
	data$Region = factor(data$Region, levels = unique(data$Region))
	# hilight
	# 1. rows, color, ... or
	# 2. chr, start, end, color
	ggs_hilight = list()
	hinames = names(hilights)
	if (!is.null(hinames)) {
		hilights = list(hilights)
	}
	hdata = NULL
	for (hilight in hilights) {
		if ('snps' %in% names(hilight)) {
			hdata_tmp = data[Snp %in% hilight$snps]
			hdata = ifelse(is.null(hdata), hdata_tmp, rbind(hdata, hdata_tmp))
			rm(hdata_tmp)
		} else if ('chr' %in% names(hilight)) {
			hdata_tmp = data[Chr == hilight$chr & Pos >= as.numeric(hilight$start) & Pos <= as.numeric(hilight$end)]
			hdata = ifelse(is.null(hdata), hdata_tmp, rbind(hdata, hdata_tmp))
			rm(hdata_tmp)
		} else {
			stop('Either "snps" or "chr" needed to locate the snps to be highlighted.')
		}
	}
	if (!is.null(hdata)) {
		hicolor = ifelse(is.null(hilight$color), "red", hilight$color)
		ggs_hilabel = list()
		if (hilabel) {
			ggs_hilabel = list(geom_label_repel = list(aes_string(x = 'X', y = 'Y', label = 'Snp'), color = "white", fill = hicolor, data = hdata, inherit.aes = FALSE, alpha = .8))
		}
		ggs_hilight = c(
			ggs_hilight,
			ggs_hilabel,
			list(geom_point = list(aes_string(x = 'X', y = 'Y'), color = hicolor, data = hdata, inherit.aes = FALSE))
		)
	}
	ggs = c(list(
		expand_limits      = list(y = c(0, max(data[, 'Y']) * 1.05)),
		geom_point         = list(aes_string(color = 'Region'), alpha = .8),
		guides             = list(color = FALSE),
		scale_color_manual = list(values = rep(c("grey", "skyblue"), nrow(rdata))),
		scale_y_continuous = list(expand = c(0, 0)),
		scale_x_continuous = list(expand = expand_scale(mult = c(0.01, 0.01)), label = rdata$Region, breaks= rdata$centerPos),
		theme_bw           = list(),
		theme              = list(
			panel.grid.minor   = element_blank(),
			panel.grid.major.x = element_blank(),
			axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5)
		),
		labs               = list(x = "", y = "-log10(p-value)")
		), ggs_hilight, ggs)
	plot.xy (data, plotfile, x = 'X', y = 'Y', ggs = ggs, devpars = devpars)
}

plot.pairs = function(data, plotfile = NULL, params = list(), ggs = list(),
	devpars = list(res = 300, width = 400, height = 400)) {

	library(GGally)
	p = do.call(ggpairs, c(list(data), params))
	for (n in names(ggs)) {
		p = p + do.call(n, ggs[[n]])
	}
	ncols          = ifelse(is.null(params$columns), ncol(data), length(params$columns))
	devpars$width  = devpars$width*ncols
	devpars$height = devpars$height*ncols
	save.plot(p, plotfile, devpars)
}

.surv.cut = function(data, vars = NULL, by = 'maxstat', labels = c("low", "high"), keep.raw = FALSE) {
	# by = mean, median, q25, q75, asis, maxstat
	# or a value as cutoff
	cnames = colnames(data)
	if (is.null(vars)) {
		vars = c(cnames[3:ncol(data)])
	}
	p = list()
	for (v in vars) {
		vdata = as.vector(unlist(data[,v]))
		if (keep.raw) {
			rawv  = paste0('raw.', v)
			data[, rawv] = data[, v]
		}
		if (by == 'asis' || !is.numeric(vdata) || length(levels(factor(vdata))) <= 5) {
			data[, v] = as.factor(as.character(vdata))
		} else if (by == 'maxstat') {
			data.cut = surv_cutpoint(data, time = cnames[1], event = cnames[2], variables = v)
			data.cat = surv_categorize(data.cut)
			data[, v] = data.cat[, v]
			p[[v]] = plot(data.cut)[[v]]
		} else {
			if (by == 'mean') {
				vcut = mean(vdata)
			} else if (by == 'median') {
				vcut = median(vdata)
			} else if (by == 'q25') {
				vcut = quantile(vdata, .25)
			} else if (by == 'q75') {
				vcut = quantile(vdata, .75)
			} else {
				vcut = by
			}
			icut = as.integer(median(which(vdata %in% vcut)))
			tmp = sapply(1:length(vdata), function(i) {
				if (vdata[i] < vcut) {
					labels[1]
				} else if (vdata[i] > vcut) {
					labels[2]
				} else if (i < icut) {
					labels[1]
				} else {
					labels[2]
				}
			})
			data[,v] = tmp
		}
	}
	# plot is a list of plots named with the variables
	return(list(data = data, plot = p))
}

plot.survival.km = function(data, plotfile = NULL, params = list(
	risk.table  = TRUE,
	pval        = TRUE,
	cut         = 'maxstat',
	cutlabels   = c('low', 'high')
), ggs = list(table = NULL), devpars = DEVPARS) {
	# @description:
	#   Use Kaplan Merier to do survival analysis
	#   See: http://www.sthda.com/english/wiki/survival-analysis-basics#kaplan-meier-survival-estimate
	# @params:
	#   dat: The data frame, where:
	#     - 1st col is time
	#     - 2nd col is status
	#     - 3rd col is the variable (have to be categorical)
	#   varname: The variable name to display instead of var (3rd colname)
	#   params:  The params used for ggsurvplot. If provided, global plot$params will not be override
	# @returns:
	#   list(plot, table, cutplot)
	library(survival)  # 2.44-1.1
	library(survminer) # 0.4.4.999

	ret = list()
	cnames  = colnames(data)
	variable = cnames[3]

	cutmethod = list.get(params, 'cut', 'maxstat')
	cutlabels = list.get(params, 'cutlabels', c('low', 'high'))

	params$cut       = NULL
	params$cutlabels = NULL

	survcut = .surv.cut(data, by = cutmethod, labels = cutlabels, keep.raw = TRUE)
	data    = survcut$data

	ret$cutplot = ifelse(length(survcut$plot) > 0, survcut$plot[[variable]], NULL)
	ret$data    = data

	# compose formula
	fmula  = sprintf('Surv(%s, %s) ~ %s',
		bQuote(cnames[1]),
		bQuote(cnames[2]),
		bQuote(variable))
	modfit = surv_fit(as.formula(fmula), data = data)

	returnTable = list.get(params, 'returnTable', FALSE)
	params$returnTable = NULL
	params$risk.table = list.get(params, 'risk.table', TRUE)
	params$pval = list.get(params, 'pval', TRUE)

	ret$plot = do.call(ggsurvplot, c(list(modfit, data = data), params))
	if (length(levels(as.factor(data[, 3]))) < 2) {
		ret$table = data.frame(
			var    = variable,
			method = 'Kaplan Merier',
			test   = 'LogRank',
			groups = 1,
			df     = 1,
			pvalue = '1.00E+00'
		)
	} else {
		counts = as.data.frame(table(data[, variable]))
		ret$table = data.frame(
			var    = variable,
			method = 'Kaplan Merier',
			test   = 'LogRank',
			groups = paste(apply(counts, 1, function(row) paste(row, collapse = ':')),
				collapse = ', '),
			df     = 1,
			pvalue = sprintf('%.2E', surv_pvalue(modfit, data = data)$pval)
		)
	}

	tableggs       = list.get(ggs, 'table', list())
	ggs$table      = NULL
	ret$plot$table = apply.ggs(ret$plot$table, tableggs)
	ret$plot$plot  = apply.ggs(ret$plot$plot, ggs)
	save.plot(ret$plot, plotfile, devpars)
	if (!is.null(plotfile) && !is.null(ret$cutplot)) {
		save.plot(
			ret$cutplot,
			sprintf("%s.cut.png", tools::file_path_sans_ext(plotfile)),
			devpars)
	}
	return(ret)
}


plot.survival.cox = function(data, plotfile = NULL, params = list(
	risk.table  = TRUE,
	pval        = TRUE,
	# mean, median, q25, q75, asis, or a value as cutoff
	# asis uses the levels in 3rd column, requring categorical values
	# or if 3rd column has <= 5 levels, regard it as categorical,
	# ignore this argument
	cut = 'maxstat',
	cutlabels = c('low', 'high')
), ggs = list(table = NULL), devpars = DEVPARS) {
	# @description:
	#   Use Kaplan Merier to do survival analysis
	#   See: http://www.sthda.com/english/wiki/cox-proportional-hazards-model
	# @params:
	#   dat: The data frame, where:
	#     - 1st col is time
	#     - 2nd col is status (censoring)
	#     - 3rd col is the target variable
	#     - ... covariates ...
	#   varname: The variable name to display instead of var (3rd colname)
	#   params:  The params used for ggsurvplot. If provided, global plot$params will not be override
	# @returns:
	#   list(data, plot, cutplot, forest, table)
	library(survival)
	library(survminer)

	ret      = list()
	cnames   = colnames(data)
	variable = cnames[3]
	if (length(cnames) > 3) {
		for (i in 4:length(cnames)) {
			data[, i] = as.numeric(as.factor(data[, i]))
		}
	}

	# compose formula
	fmula    = sprintf('Surv(%s, %s) ~ .', bQuote(cnames[1]), bQuote(cnames[2]))
	# do cox regression
	fit0     = coxph(as.formula(fmula), data = data)
	sumfit    = summary(fit0)
	ret$table = data.frame(
		var    = rownames(sumfit$conf.int),
		method = 'Cox Regression',
		test   = 'Log Rank')

	ret$table         = cbind(
		ret$table, sumfit$conf.int[, -2, drop = FALSE], sumfit$coefficients[, 4:5, drop = FALSE])
	ret$table$df      = sumfit$logtest[2]
	ret$table$modpval = sumfit$logtest[3]
	ret$forest        = ggforest(fit0, data = data)
	colnames(ret$table)[4:8] = c('HR', 'CI95_1', 'CI95_2', 'z', 'pvalue')

	# now we cut the variable to plot it
	cutmethod = list.get(params, 'cut', 'maxstat')
	cutlabels = list.get(params, 'cutlabels', c('low', 'high'))
	params$cut = NULL
	params$cutlabels = NULL
	survcuts  = .surv.cut(data, by = cutmethod, labels = cutlabels, vars = variable, keep.raw = TRUE)
	ret$cutplot = ifelse(length(survcuts$plot) > 0, survcuts$plot[[variable]], NULL)
	ret$data    = survcuts$data
	survcuts$data[[paste0('raw.', variable)]] = NULL
	fit     = coxph(as.formula(fmula), data = survcuts$data)
	groups  = levels(as.factor(survcuts$data[,variable]))
	newdata = data.frame(groups = groups)
	colnames(newdata) = variable
	if (length(cnames) > 3) {
		for (i in 4:length(cnames)) {
			var = cnames[i]
			newdata[[var]] = mean(as.numeric(survcuts$data[,var]), na.rm = TRUE)
		}
	}
	fit2 = survfit(fit, newdata = newdata)
	risktable = list.get(params, 'risk.table', TRUE)
	params$risk.table = FALSE
	pval = list.get(params, 'pval', TRUE)
	if (pval) {
		params$pval = sprintf(
			'P = %.3E\nHR = %.3f (95%% CI %.3f ~ %.3f)',
			sumfit$logtest[3], ret$table[1,4], ret$table[1,5], ret$table[1,6])
	}

	ret$plot = do.call(ggsurvplot,
		c(list(fit2, data = newdata, legend.labs=paste(variable, '=', groups)), params))

	# need to hack a risk table.
	# See https://github.com/kassambara/survminer/issues/231#issuecomment-309453646
	if (risktable) {
		fmula = sprintf('Surv(%s, %s) ~ %s',
				bQuote(cnames[1]), bQuote(cnames[2]), bQuote(variable))
		kmfit = surv_fit(as.formula(fmula), data = survcuts$data)
		rtparams = list(kmfit, risk.table = TRUE, data = survcuts$data)
		if (!is.null(params$break.time.by)) {
			rtparams = c(rtparams, list(break.time.by = params$break.time.by))
		}
		if (!is.null(params$xlim)) {
			rtparams = c(rtparams, list(xlim = params$xlim))
		}
		kmplot = do.call(ggsurvplot, rtparams)
		kmplot$plot = ret$plot$plot
		ret$plot = kmplot
		tableggs = list.get(ggs, 'table', list())
		ggs$table = NULL
		ret$plot$table = apply.ggs(ret$plot$table, tableggs)
	}
	ret$plot$plot = apply.ggs(ret$plot$plot, ggs)
	save.plot(ret$plot, plotfile, devpars)
	if (!is.null(plotfile) && !is.null(ret$cutplot)) {
		save.plot(
			ret$cutplot[[variable]],
			sprintf("%s.cut.png", tools::file_path_sans_ext(plotfile)),
			devpars)
	}
	if (!is.null(plotfile)) {
		save.plot(
			ret$forest,
			sprintf("%s.forest.png", tools::file_path_sans_ext(plotfile)),
			devpars)
	}
	return(ret)
}

# The data should be a data.frame with rows as samples and columns as genes
#         Group   Gene1  Gene2 ... GeneN
# sample1 Control 1      2     ... x
# ...
# params$has.group indicates if the data has first column as group or not
plot.mds = function(data, plotfile, params = list(has.group = TRUE, ndim = 2),
	ggs = list(), devpars = DEVPARS) {

	data = as.data.frame(data)
	groupdata = NULL
	if (list.get(params, 'has.group', TRUE)) {
		groupdata = data[, 1, drop = FALSE]
		data = data[, -1, drop = FALSE]
	}
	plotdata = as.data.frame(cmdscale(dist(scale(data)), k = list.get(params, 'ndim', 2)))
	if(!is.null(groupdata)) {
		plotdata = cbind(plotdata,Group = groupdata[rownames(plotdata),])
		hulldata = NULL
		for (grup in levels(as.factor(plotdata$Group))) {
			if (is.null(hulldata)) {
				hulldata = plotdata[plotdata$Group == grup,,drop = FALSE][
					chull(plotdata[plotdata$Group == grup, 1:2]), ]
			} else {
				hulldata = rbind(hulldata, plotdata[plotdata$Group == grup,,drop = FALSE][
					chull(plotdata[plotdata$Group == grup, 1:2]), ])
			}
		}
		print(hulldata)
		ggs = c(list(
			geom_point = list(aes_string(color = 'Group'), shape = 21),
			geom_polygon = list(aes_string(fill = 'Group', group = 'Group'), alpha = .3, data = hulldata),
			geom_text_repel = list(label = rownames(plotdata), min.segment.length = unit(0, 'lines')),
			theme_bw = list(),
			labs = list(x = 'Dimension 1', y = 'Dimension 2'),
			theme = list(axis.text.x = element_blank(),  # remove x-axis text
				axis.text.y = element_blank(), # remove y-axis text
				axis.ticks = element_blank(),  # remove axis ticks
				panel.background = element_blank(),
				panel.grid.major = element_blank(),  #remove major-grid labels
				panel.grid.minor = element_blank(),  #remove minor-grid labels
				plot.background = element_blank())
		), ggs)
	} else {
		ggs = c(list(
			geom_point = list(shape = 21),
			geom_text_repel = list(label = rownames(plotdata), min.segment.length = unit(0, 'lines')),
			theme_bw = list(),
			labs = list(x = 'Dimension 1', y = 'Dimension 2'),
			theme = list(axis.text.x = element_blank(),  # remove x-axis text
				axis.text.y = element_blank(), # remove y-axis text
				axis.ticks = element_blank(),  # remove axis ticks
				panel.background = element_blank(),
				panel.grid.major = element_blank(),  #remove major-grid labels
				panel.grid.minor = element_blank(),  #remove minor-grid labels
				plot.background = element_blank())
		), ggs)
	}
	plot.xy(plotdata, plotfile, x = 'V1', y = 'V2', ggs = ggs, devpars = devpars)

}

