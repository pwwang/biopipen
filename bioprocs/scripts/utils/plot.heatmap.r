require('ggplot2')
require('ggdendro')
require('gtable')
require('grid')
if (!exists('plotHeatmap')) {
	# dendro: 
	#   TRUE/'both': drow dendro for both dimension
	#   'row': only draw dendro for rows
	#   'col': only draw dendro for cols
	# rows:
	#   Order of rows (rownames) if dendro is not drawn for rows
	# cols:
	#   Order of cols (colnames) if dendro is not drawn for cols
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
}