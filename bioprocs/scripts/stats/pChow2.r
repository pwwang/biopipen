library(methods)
{{rimport}}('plot.r', '__init__.r')

infile   = {{i.infile | R}}
outfile  = {{o.outfile | R}}
outdir   = {{o.outdir | R}}
plotchow = {{args.plot | R}}
devpars  = {{args.devpars | R}}
fdr      = {{args.fdr | R}}
inopts   = {{args.inopts | R}}
covfile  = {{args.cov | R}}

if (covfile != '') {
	cov    = read.table(covfile, header = T, row.names = 1, sep = "\t", check.names = F)
	#indata = cbind(cov, indata)
}

onecase = function(indata, case) {
	ncols.nocov = ncol(indata)

	subgroups    = levels(factor(indata[, ncols.nocov]))
	nesamples    = c() # groups with not enough samples
	for (g in subgroups) {
		gdata = indata[which(indata[, ncols.nocov] == g),,drop = F]
		if (nrow(gdata) < 2)
			nesamples = c(nesamples, rownames(gdata))
	}
	if (length(nesamples) > 0) {
		indata = indata[setdiff(rownames(indata), nesamples),,drop = F]
		rm(nesamples)
	}
	subgroups    = levels(factor(indata[, ncols.nocov]))

	ngroups      = length(subgroups)
	cnames.nocov = colnames(indata)
	# for plotting
	Y            = cnames.nocov[ncols.nocov - 1]
	X            = cnames.nocov[ncols.nocov - 2]
	Group        = cnames.nocov[ncols.nocov]

	if (covfile != '') {
		indata = cbind.fill(cov, indata)
	}

	cnames = colnames(indata)
	ncols  = ncol(indata)
	fmula  = as.formula(sprintf('%s ~ %s', Y, paste(bQuote(cnames[1:(ncols-2)]), collapse = '+')))
	# do regression on pooled:
	lm_p  = lm(fmula, data = indata)
	ssr_p = sum(lm_p$residuals ^ 2)

	ssr_groups = 0
	# do for subgroups
	lms = list()
	lmfile = paste0(tools::file_path_sans_ext(outfile), '.', case, '.lms')
	lmsum  = NULL
	for (i in 1:ngroups) {
		m = lm(fmula, data = indata[which(indata[, ncols] == subgroups[i]), , drop = F])
		ssr_groups = ssr_groups + sum(m$residuals ^ 2)
		lms[[subgroups[i]]] = m
		m_sum = summary(m)
		m_sum$coefficients = cbind(m_sum$coefficients, Group = subgroups[i])
		lmsum = rbind(lmsum, m_sum$coefficients)
	}
	write.table(pretty.numbers(lmsum, list(
		`Estimate..Std. Error..t value` = '%.3f',
		`Pr(>|t|)` = '%.3E'
	)), lmfile, col.names = T, row.names = T, sep = "\t", quote = F)
	K      = ncols.nocov - 1
	J      = (ngroups - 1) * K
	DF     = nrow(indata) - ngroups * K
	fscore = (ssr_p - ssr_groups) * DF / ssr_groups / J
	pvalue = pf(fscore, J, DF, lower.tail = F)

	# save the results
	out = data.frame(Case = case, Groups = ngroups, K = K, J = J, DF = DF, fstat = fscore, pval = pvalue)

	if (!is.na(pvalue) && (plotchow == T || pvalue < plotchow)) {
		# pooled regression line
		model2eq = function(model) {
			vars = colnames(model$model)
			cf   = round(model$coefficients, 2)
			paste0(vars[1], ' = ', cf[1], paste(
				sapply(
					2:length(cf),
					function(x) paste0(
						if(is.na(cf[x]))''else if(cf[x]>=0)'+'else'-', abs(cf[x]), '*', vars[x])
				), collapse = ''
			))
		}
		labels = c()
		for (subgroup in subgroups) {
			m = lms[[subgroup]]
			labels = c(labels, paste0(subgroup, ': ', model2eq(m)))
		}
		# pooled
		labels = c(labels, paste0('Pooled: ', model2eq(lm_p)))
		ggs = list(
			geom_smooth = list(
				aes(color = 'Pooled'),
				method   = 'lm',
				se       = F,
				linetype = "twodash"
			),
			geom_smooth = list(
				aes_string(color = Group),
				method  = 'lm',
				se      = F
			),
			theme = list(
				legend.position = "bottom"
			),
			guides = list(
				color=guide_legend(ncol=1),
				shape=F
			),
			scale_color_manual = list(
				values = c(scales::hue_pal()(ngroups), '#555555'), 
				name   = Group,
				limit  = c(subgroups, 'Pooled'),
				labels = c(labels, 'Pooled')
			)
		)
		plotfile = paste0(tools::file_path_sans_ext(outfile), '.', case, '.png')
		plot.scatter(
			indata, 
			plotfile, 
			x      = ncols-2,
			y      = ncols-1,
			ggs    = ggs,
			params = list(aes_string(shape = Group, color = Group))
		)
	}
	return(out)
}

indata  = read.table.inopts(infile, inopts)
cases   = levels(factor(indata$Case))
incols  = colnames(indata)
results = NULL
for (case in cases) {
	log2pyppl('- doing case ', case, '...')
	rmcol = c('Case')
	indata2 = indata[which(indata$Case == case),, drop = F]
	if ('Sample' %in% incols) {
		rownames(indata2) = indata2$Sample
		rmcol = c(rmcol, 'Sample')
	}
	cols = setdiff(incols, rmcol)
	out  = onecase(indata2[,  cols, drop = F], case)
	results = if (is.null(results)) out else rbind(results, out)
}
if (fdr != F) {
	if (fdr == T) fdr = 'BH'
	results$fdr = p.adjust(results$pval, method = fdr)
}

write.table(pretty.numbers(results, list(
	fstat = '%.3f',
	pval  = '%.3E',
	fdr   = '%.3E'
)), outfile, col.names = T, row.names = F, sep = "\t", quote = F)
	