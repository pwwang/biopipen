require('ggplot2')
require('ggdendro')
require('gtable')
require('grid')
if (!exists('plotHeatmap')) {
	plotHeatmap = function (m, filename, ggs = list(), devpars = list(res=300, width=2000, height=2000), dendro = TRUE) {
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
		
		x = as.matrix(m)
		# Dendrogram 1
		ddrleft   = hclust(dist(x))
		plotdleft = ggplot(segment(dendro_data(as.dendrogram(ddrleft)))) + 
		geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
		theme_none + theme(plot.margin = margin(r=-5, l = 10)) + 
		coord_flip() + 
		scale_y_reverse(expand = c(0, 0)) +
		scale_x_discrete(expand = c(0, 0.5)) 
		
		# Dendrogram 2
		ddrtop   = hclust(dist(t(x)))
		plotdtop = ggplot(segment(dendro_data(as.dendrogram(ddrtop)))) + 
		geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
		theme_none + theme(plot.margin = margin(b=-15, t = 10)) +
		scale_x_discrete(expand = c(0, 0.5)) +
		scale_y_continuous(expand = c(0, 0))
		
		
		mm   = stack(as.data.frame(m))
		rows = rownames(m)
		if (is.null(rows)) {
		rows = paste('ROW', ddrleft$order)
		}
		#colidx  = match(1:ncol(m), ddrtop$labels$label)
		colidx  = ddrtop$order
		rowidx  = ddrleft$order
		
		mm$rows = rows
		#[as.numeric(levels(ddrleft$labels$label)), as.numeric(levels(ddrtop$labels$label))]
		plothm  = ggplot(mm, aes(ind, rows)) + geom_tile(aes(fill=values)) +
		xlim(unique(as.vector(mm$ind))[colidx]) + scale_y_discrete(position="right", limits=mm$rows[rowidx])
		
		for (i in 1:length(ggplus)) {
		plothm = plothm + ggplus[[i]]
		}
		if (length(ggs) > 0) {
		for (i in 1:length(ggs)) {
			plothm = plothm + ggs[[i]]
		}
		}
		
		if (dendro) {
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
		} else {
			print(plothm)
		}
		dev.off()
	}
}