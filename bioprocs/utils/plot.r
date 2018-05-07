require('ggplot2')
require('ggdendro')
require('gtable')
require('grid')
require('ggpubr')

# see http://www.sthda.com/english/rpkgs/ggpubr/reference/ggscatter.html for all parameters
scatter = function(data, plotfile, x = 1, y = 2, params = list(), devpars = list(res = 300, width = 2000, height = 2000)) {
	do.call(png, c(list(filename=plotfile), devpars))
	cnames = colnames(data)
	if (is.integer(x)) x = cnames[x]
	if (is.integer(y)) y = cnames[y]
	params$data = data
	params$x    = x
	params$y    = y
	# put cor.coeff label in the center if cor.coeff.args$label.x is not given
	if ('cor.coeff.args' %in% names(params) && !'label.x' %in% names(params$cor.coeff.args))
		params$cor.coeff.args$label.x = (max(data[,x]) + min(data[,x])) / 2
	print(do.call(ggscatter, params))
	dev.off()
}


plotBoxplot = function (m, filename, ggs = list(), devpars = list(res=300, width=2000, height=2000)) {
	do.call(png, c(list(filename=filename), devpars))

	ggplus = list(
		ylab("Value"),
		theme(axis.title.x = element_blank()),
		theme(axis.text.x  = element_text(angle = 60, hjust = 1)),
		theme(plot.margin  = unit(c(1,1,1,1), "cm"))
	)
	plotout = ggplot(stack(as.data.frame(m)), aes(ind, values)) + geom_boxplot()
	for (i in 1:length(ggplus)) {
		plotout = plotout + ggplus[[i]]
	}
	if (length(ggs) > 0) {
		for (i in 1:length(ggs)) {
			plotout = plotout + ggs[[i]]
		}
	}
	print(plotout)
	dev.off()
}

plotHeatmap = function (m, filename, ggs = list(), devpars = list(res=300, width=2000, height=2000), dendro = TRUE, rows = NULL, cols = NULL) {
	do.call(png, c(list(filename=filename), devpars))

	theme_none <- theme(
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
		theme(axis.ticks = element_blank()),
		scale_fill_gradient2(),
		theme(axis.title.x = element_blank()),
		theme(axis.title.y = element_blank()),
		theme(axis.text.x  = element_text(angle = -60, hjust = 0)),
		#scale_y_discrete(position="right"),
		theme(legend.title = element_blank())
	)

	x      = as.matrix(m)
	rnames = rownames(m)
	cnames = colnames(m)
	if (is.null(rnames)) rnames = paste('ROW', 1:nrow(m), sep='')
	if (is.null(cnames)) cnames = paste('COL', 1:ncol(m), sep='')
	rownames(m) = rnames
	colnames(m) = cnames
	# Dendrogram 1
	if (dendro == TRUE || dendro == 'both' || dendro == 'row') {
		ddrleft   = hclust(dist(x))
		plotdleft = ggplot(segment(dendro_data(as.dendrogram(ddrleft)))) +
		geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
		theme_none + theme(plot.margin = margin(r=-5, l = 10)) +
		coord_flip() +
		scale_y_reverse(expand = c(0, 0)) +
		scale_x_discrete(expand = c(0, 0.5))
		rowidx  = ddrleft$order
	} else {
		plotdleft = NULL
		rowidx    = if (!is.null(rows)) match(rows, rnames) else 1:length(rnames)
		rowidx    = rev(rowidx)
	}

	# Dendrogram 2
	if (dendro == TRUE || dendro == 'both' || dendro == 'col') {
		ddrtop   = hclust(dist(t(x)))
		plotdtop = ggplot(segment(dendro_data(as.dendrogram(ddrtop)))) +
		geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
		theme_none + theme(plot.margin = margin(b=-15, t = 10)) +
		scale_x_discrete(expand = c(0, 0.5)) +
		scale_y_continuous(expand = c(0, 0))
		colidx  = ddrtop$order
	} else {
		plotdtop = NULL
		colidx   = if (!is.null(cols)) match(cols, cnames) else 1:length(cnames)
	}

	mm        = stack(as.data.frame(m[rowidx, colidx]))
	mm$rnames = rnames
	plothm    = ggplot(mm, aes(ind, rnames)) + geom_tile(aes(fill=values)) +
	xlim(unique(as.vector(mm$ind))) + scale_y_discrete(position="right", limits=rnames)

	for (i in 1:length(ggplus)) {
		plothm = plothm + ggplus[[i]]
	}
	if (length(ggs) > 0) {
		for (i in 1:length(ggs)) {
			plothm = plothm + ggs[[i]]
		}
	}

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

plotHist = function (m, filename, ggs = list(), devpars = list(res=300, width=2000, height=2000)) {
	do.call(png, c(list(filename=filename), devpars))

	ggplus  = list(
		xlab("Value"),
		ylab("Count")
	)
	plotout = ggplot(stack(as.data.frame(t(m)))) + geom_histogram(aes(values))
	for (i in 1:length(ggplus)) {
		plotout = plotout + ggplus[[i]]
	}
	if (length(ggs) > 0) {
		for (i in 1:length(ggs)) {
			plotout = plotout + ggs[[i]]
		}
	}
	print(plotout)
	dev.off()
}

plotMAplot = function (m, filename, ggs = list(), devpars = list(res=300, width=2000, height=2000)) {
	do.call(png, c(list(filename=filename), devpars))
	ggplus = list(
		theme(legend.position = "none"),
		geom_hline(color = "blue3", yintercept=0, linetype="dashed"),
		stat_smooth(se = FALSE, method = "loess", color = "red3")
	)

	cnames = names(m)
	A      = if ("A" %in% cnames) m$A else m[, 1]
	M      = if ("M" %in% cnames) m$M else m[, 2]
	thres  = if ("threshold" %in% cnames) m$threshold else if (ncol(m) > 2) m[, 3] else NULL

	if (is.null(thres)) {
		plotout = ggplot(m, aes(x = A, y = M)) + geom_point(size = 1.5, alpha = 1/5)
	} else {
		plotout = ggplot(m, aes(x = A, y = M)) + geom_point(size = 1.5, alpha = 1/5, aes(color=factor(thres)))
	}

	for (i in 1:length(ggplus)) {
		plotout = plotout + ggplus[[i]]
	}
	if (length(ggs) > 0) {
		for (i in 1:length(ggs)) {
			plotout = plotout + ggs[[i]]
		}
	}
	print(plotout)
	dev.off()
}

plotPie = function (m, filename, ggs = list(), devpars = list(res=300, width=2000, height=2000)) {
	percent = function(x) paste0(format(round(x*100, 1), nsmall = 1), "%")
	do.call(png, c(list(filename=filename), devpars))
	N         = nrow(m)
	Group     = colnames(m)
	oValue    = if (N == 1) m[1,] else colSums(m)
	oValue    = as.vector(as.matrix(oValue))
	GroupName = paste0(Group, " (", oValue, ")")
	Value     = 100*oValue/sum(oValue)
	Labels    = percent(Value / 100)
	df        = data.frame (Group = Group, Value = Value)
	df$Group  = factor(df$Group, levels = df$Group)
	#if (length(Group) == 2) {
	#	ytext = 100 - cumsum(rev(Value)) + rev(Value) / 2
	#} else {
		ytext = cumsum(rev(Value)) - rev(Value) / 2
	#}
	pie =   ggplot(df, aes(x="", y=Value, fill=Group)) +
			geom_bar(width = 1, stat = "identity") +
			geom_text(aes(y = ytext), label = rev(Labels)) +
			coord_polar("y") +
			scale_fill_discrete(labels = GroupName, breaks = Group)

	if (length(ggs) > 0) {
		for (i in 1:length(ggs)) {
			pie = pie + ggs[[i]]
		}
	}
	print(pie)
	dev.off()
}

plotUpset = function(mat, filename, params, devpars) {
	library(UpSetR)
	do.call(png, c(list(filename = filename), devpars))
	png (filename, res=300, width=2000, height=2000)
	if (! "nintersects" %in% names(params)) {
		params$nintersects = NA
	}
	v = do.call (upset, c(list(data = mat, sets = colnames(mat)), params))
	print (v)
	dev.off()
}

plotVenn = function(mat, filename, params, devpars) {
	library(VennDiagram)
	rnames = rownames(mat)
	if (is.null(rnames)) {
		rnames = paste('ROW', 1:nrow(mat), sep = '')
		rownames(mat) = rnames
	}
	x = list()
	for (cname in colnames(mat)) {
		x[[cname]] = rownames(mat[which(mat[, cname] == 1), cname, drop=F])
	}
	ps = list(x = x, filename = filename, height = devpars$height, width = devpars$width, resolution = devpars$res, imagetype = "png", alpha = .5, fill = rainbow(ncol(mat)))

	do.call(venn.diagram, c(ps, params))

}

plotVolplot = function(m, filename, ggs = list(), devpars = list(res=300, width=2000, height=2000)) {
	do.call(png, c(list(filename=filename), devpars))
	m = as.data.frame(m)

	cnames = names(m)
	logfc  = if ("logFC" %in% cnames) m$logFC else m[, 1]
	fdr    = if ("FDR" %in% cnames) m$FDR else m[, 2]
	fdr    = -log10(fdr)

	# cutoffs
	logfccut    = if ("logFCCut" %in% cnames) m$logFCCut[1] else 2
	fdrcut      = if ("FDRCut" %in% cnames) m$FDRCut[1] else 0.05
	fdrcutlabel = fdrcut
	fdrcut      = -log10(fdrcut)

	threshold = as.factor(abs(logfc) > logfccut & fdr > fdrcut)

	xm = min(max(abs(logfc)), 10)
	ym = min(max(fdr), 15)
	if (xm <= logfccut) logfccut = 1

	df  = data.frame(logfc, fdr)
	plotout = ggplot(data=df, aes(x=logfc, y=fdr)) +
	  geom_point(alpha=0.4, size=1.75, aes(color=threshold)) +
	  geom_hline(yintercept=fdrcut, linetype="dashed", color="blue3") +
	  geom_vline(xintercept=c(-logfccut, logfccut), linetype="dashed", color = "red3") +
	  xlim(c(-xm, xm)) + ylim(c(0, ym)) +
	  geom_text(aes(xm, fdrcut, label = paste("p", "=", fdrcutlabel), vjust = -1, hjust = 1), color="blue3") +
	  geom_text(aes(-logfccut, ym, label = paste(-logfccut, "fold"),  vjust = 1, hjust = -0.1), color="red3") +
	  geom_text(aes(+logfccut, ym, label = paste(paste('+', logfccut, sep=""), "fold"),  vjust = 1, hjust = -0.1), color="red3")

	ggplus = list(
		theme(legend.position = "none"),
		xlab("log2 Fold Change"),
		ylab("-log10(Pval)")
	)

	for (i in 1:length(ggplus)) {
		plotout = plotout + ggplus[[i]]
	}
	if (length(ggs) > 0) {
		for (i in 1:length(ggs)) {
			plotout = plotout + ggs[[i]]
		}
	}
	print(plotout)
	dev.off()
}
