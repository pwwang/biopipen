library(plotROC)

cnames = as.logical({{args.cnames | R}})
rnames = as.logical({{args.rnames | R}})
df = read.table({{in.infile | quote}}, sep = "\t", check.names = F, header = cnames, row.names = if (rnames) 1 else NULL)

if (!cnames) {
	colnames(df) = c('D', paste0('M', 1:(ncol(df)-1)))
}
headers = colnames(df)

if (ncol(df) == 2) {
	rocfile = file.path({{out.outdir | quote}}, 'roc.png')
	aucfile = file.path({{out.outdir | quote}}, 'auc.txt')
	baseroc = ggplot(df, aes(d = df[, 1], m = df[, 2]))
	baseroc = baseroc + do.call(geom_roc, {{args.params | Rlist}})
	auc     = matrix(round(calc_auc(baseroc)$AUC, 3), ncol=1)

	{% if args.showauc %}
	baseroc = baseroc + annotate("text", x = .9, y = .1, label = paste('AUC', auc, sep=': '))
	{% endif %}

	ggs     = {{args.ggs | Rlist}}
	for (i in 1:length(ggs)) {
		baseroc = baseroc + ggs[[i]]
	}

	do.call(png, c(filename = rocfile, {{args.devpars | Rlist}}))
	print(baseroc)
	dev.off()

	rownames(auc) = headers[-1]
	colnames(auc) = 'AUC'
	write.table(auc, aucfile, sep = "\t", quote = F)
	
} else if (ncol(df) > 2) {
{% 	if args.combine %}
	rocfile = file.path({{out.outdir | quote}}, 'roc.png')
	aucfile = file.path({{out.outdir | quote}}, 'auc.txt')
	rocdf   = melt_roc(df, 1, 2:ncol(df))
	baseroc = ggplot(rocdf, aes(d = D, m = M, color = name))
	baseroc = baseroc + do.call(geom_roc, {{args.params | Rlist}})
	auc     = matrix(round(calc_auc(baseroc)$AUC, 3), ncol = 1)

	{% if args.showauc %}
	auclabels = NULL
	for (i in 1:length(headers[-1])) {
		auclabels = c(auclabels, paste(headers[-1][i], auc[i, 1], sep = ': '))
	}
	baseroc = baseroc + scale_color_discrete(name = "AUC", breaks = headers[-1], labels = auclabels)
	{% 	endif %}

	ggs     = {{args.ggs | Rlist}}
	for (i in 1:length(ggs)) {
		baseroc = baseroc + ggs[[i]]
	}
	do.call(png, c(filename = rocfile, {{args.devpars | Rlist}}))
	print(baseroc)
	dev.off()

	rownames(auc) = headers[-1]
	colnames(auc) = 'AUC'
	write.table(auc, aucfile, sep = "\t", quote = F)

{% 	else %}
	aucfile = file.path({{out.outdir | quote}}, 'auc.txt')
	aucs    = matrix(NA, ncol = 1, nrow = length(headers) - 1)
	rownames(aucs)   = headers[-1]
	colnames(aucs)   = 'AUC'
	for (m in headers[-1]) {
		rocfile = file.path({{out.outdir | quote}}, paste0(m, '.roc.png'))
		rocdf = df[, c(headers[1], m)]
		baseroc = ggplot(rocdf, aes(d = df[, 1], m = df[, m]))
		baseroc = baseroc + do.call(geom_roc, {{args.params | Rlist}})
		auc     = matrix(round(calc_auc(baseroc)$AUC, 3), ncol=1)
		aucs[m, 1] = auc
		{% if args.showauc %}
		baseroc = baseroc + annotate("text", x = .9, y = .1, label = paste('AUC', auc, sep=': '))
		{% endif %}

		ggs     = {{args.ggs | Rlist}}
		for (i in 1:length(ggs)) {
			baseroc = baseroc + ggs[[i]]
		}
		do.call(png, c(filename = rocfile, {{args.devpars | Rlist}}))
		print(baseroc)
		dev.off()
	}
	write.table(aucs, aucfile, sep = "\t", quote = F)

{% 	endif %}
}