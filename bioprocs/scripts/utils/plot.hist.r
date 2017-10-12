require ('ggplot2')
if (!exists('plotHist')) {
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
}