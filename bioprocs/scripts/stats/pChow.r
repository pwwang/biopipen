library(methods)
{{rimport}}('plot.r', '__init__.r')

infile   = {{i.infile | R}}
outfile  = {{o.outfile | R}}
outdir   = {{o.outdir | R}}
plotchow = {{args.plot | R}}
devpars  = {{args.devpars | R}}
inopts   = {{args.inopts | R}}
covfile  = {{args.cov | R}}

indata = read.table(infile, header = inopts$cnames, row.names = if(inopts$rnames) 1 else NULL, sep = "\t", check.names = F)

ncols.nocov = ncol(indata)

# if you have more than 2 groups, the rest groups will do chow-test against the first one individually
subgroups    = levels(factor(indata[, ncols.nocov]))
ngroups      = length(subgroups)
cnames.nocov = colnames(indata)
# for plotting
Y            = cnames.nocov[ncols.nocov - 1]
X            = cnames.nocov[ncols.nocov - 2]
Group        = cnames.nocov[ncols.nocov]

if (covfile != '') {
	cov    = read.table(covfile, header = T, row.names = 1, sep = "\t", check.names = F)
	indata = cbind(cov, indata)
}

cnames = colnames(indata)
ncols  = ncol(indata)
fmula  = as.formula(sprintf('%s ~ %s', Y, paste(cnames[1:(ncols-2)], collapse = '+')))
# do regression on pooled:
lm_p  = lm(fmula, data = indata)
ssr_p = sum(lm_p$residuals ^ 2)

ssr_groups = 0
# do for subgroups
lms = list()
lmfile = paste0(tools::file_path_sans_ext(outfile), '.lms')
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
out = data.frame(Groups = ngroups, K = K, J = J, DF = DF, fstat = fscore, pval = pvalue)
write.table(pretty.numbers(out, list(
	fstat = '%.3f',
	pval  = '%.3E'
)), outfile, col.names = T, row.names = F, sep = "\t", quote = F)

if (plotchow == T || pvalue < plotchow) {
	# pooled regression line
	model2eq = function(model) {
		vars = colnames(model$model)
		cf   = round(model$coefficients, 2)
		paste0(vars[1], ' = ', cf[1], paste(
			sapply(
				2:length(cf),
				function(x) paste0(if(cf[x]>0)'+'else'-', cf[x], '*', vars[x])
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
	plotfile = paste0(tools::file_path_sans_ext(outfile), '.png')
	plot.scatter(
		indata, 
		plotfile, 
		x      = X,
		y      = Y,
		ggs    = ggs,
		params = list(aes_string(shape = Group, color = Group))
	)
}
