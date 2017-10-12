require ('ggplot2')
if (!exists('plotMAplot')) {
	# requires: m$A, m$M, m$threshold
	# A = mean
	# M = log2fc
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
}