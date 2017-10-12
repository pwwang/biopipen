require ('ggplot2')
if (!exists('plotBoxplot')) {
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
}