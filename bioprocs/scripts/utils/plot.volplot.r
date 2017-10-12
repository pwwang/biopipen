require('ggplot2')
if (!exists('plotVolplot')) {
	# m$logFC, m$FDR, m$logFCCut, m$FDRCut
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

}