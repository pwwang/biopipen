{{rimport}}('__init__.r', 'plot.r')

{% python from os import path %}
library(methods)
library(survival)
library(survminer)

exprfile = {{i.exprfile | quote}}
survfile = {{i.survfile | quote}}
gene     = {{i.gene | quote}}
outdir   = {{o.outdir | quote}}
lcsig    = {{args.lcsig | quote}}
inunit   = {{args.inunit | quote}}
outunit  = {{args.outunit | quote}}
params   = {{args.params | R}}
devpars  = {{args.devpars | R}}
mainggs  = {{args.ggs | R}}

prefix   = {{fn2(i.survfile), i.gene, o.outdir | lambda v1, v2, v3, path = path: path.join(v3, '{}.LCSIG{}.'.format(v1, '-' + v2 if v2 else '')) | quote }}

tableggs = mainggs$table
mainggs$table = NULL

# read raw survival data
fct = 1
{% case (args.inunit, args.outunit) %}
	{% when ('days', 'weeks') %}
	fct = 1 / 7
	{% when ('days', 'months') %}
	fct = 1 / 30
	{% when ('days', 'years') %}
	fct = 1 / 365 # 365.25?
	{% when ('weeks', 'days') %}
	fct = 7
	{% when ('weeks', 'months') %}
	fct = 7*12/365
	{% when ('weeks', 'years') %}
	fct = 7/365
	{% when ('years', 'days') %}
	fct = 365
	{% when ('years', 'weeks') %}
	fct = 365/7
	{% when ('years', 'months') %}
	fct = 12
{% endcase %}

# rows: samples
# cols: patients, time, status
survdata = read.table(survfile, header = T, row.names = 1, sep = "\t", check.names = F)
survsamples = rownames(survdata)

# read all expression data
# rows: genes
# cols: samples
exprdata = read.table.nodup(exprfile, header = T, row.names = 1, sep = "\t", check.names = F)
allgenes = rownames(exprdata)
exprsamples = colnames(exprdata)

samples = intersect(survsamples, exprsamples)
survdata = survdata[samples,, drop=F]
exprdata = exprdata[,samples, drop=F]
if (fct != 1) {
	survdata[, 1] = survdata[, 1] * fct
}

# read signature
sigenes = read.table(lcsig, header = T, row.names = 1, sep = "\t", check.names = F)
sigenes = intersect(allgenes, rownames(sigenes))

sigexpr = as.matrix(colMeans(exprdata[sigenes,, drop=F]))
colnames(sigexpr) = 'LCSIG'

if (gene != "") {
	gexpr = exprdata[gene, samples, drop=F]
	survdata = cbind(survdata[samples,], sigexpr, t(gexpr))
} else {
	survdata = cbind(survdata[samples,], sigexpr)
}

write.table(
	survdata, 
	paste0(prefix, 'rawdata.txt'),
	row.names = T,
	col.names = T,
	sep = "\t",
	quote = F
)

autogroup = function(survdata, var) {
	quantsteps = seq(0.15, 0.85, by = 0.005)
	survscores = sapply(quantsteps, function(x){
		newdata = survdata
		newdata[, var] = as.numeric(survdata[, var] >= quantile(survdata[, var], x))
		fmula = as.formula(paste('Surv(time, status) ~', var))
		tryCatch({
			coxph(fmula, data = newdata)$score
		}, error = function(x) 1)
	})
	plot.scatter(
		data.frame(step = quantsteps, survscores = survscores), 
		paste0(prefix, var, '.quant.png'), 
		ggs = list(geom_line = list()), 
		devpars = devpars
	)
	quantsteps[which.max(survscores)][1]
}

survdata[, 'LCSIG'] = c('low', 'high')[as.numeric(survdata[, 'LCSIG'] > quantile(
	survdata[, 'LCSIG'], 
	autogroup(survdata, 'LCSIG')
)) + 1]


if (gene != "") {
	survdata[, gene] = c('low', 'high')[as.numeric(survdata[, gene] > quantile(
		survdata[, gene], 
		autogroup(survdata, gene)
	)) + 1]
	Group = apply(survdata, 1, function(row) paste(row[3], row[4], sep = ':'))
	Group = matrix(Group, ncol = 1)
	colnames(Group) = paste('LCSIG', gene, sep = ':')
	survdata = cbind(survdata, Group)
}

write.table(
	survdata, 
	paste0(prefix, 'groupdata.txt'),
	row.names = T,
	col.names = T,
	sep = "\t",
	quote = F
)

plot.survival = function(data, v) {
	var            = make.names(v)
	cnames         = colnames(data)
	cnames[which(cnames == v)] = var
	colnames(data) = cnames
	varlvls        = levels(factor(data[, var]))
	params.default = list(
		legend.title = paste0('Strata(', v, ')'),
		legend.labs  = varlvls,
		xlab         = paste0("Time (", outunit ,")"),
		ylim.min     = 0.0,
		pval.coord   = c(0, 0.1)
	)
	params = update.list(params.default, params)
	fmula  = as.formula(paste('Surv(time, status) ~', var))
	modkm  = surv_fit(fmula, data = data)
	survp  = surv_pvalue(modkm)

	params$legend.title = gsub('{var}', v, params$legend.title, fixed = T)
	if (!is.null(params$pval) && params$pval != F) {
		params$pval = gsub('{pval}', sprintf('%.2E', survp$pval), params$pval, fixed = T)
	}
	ylim.min        = params$ylim.min
	params$ylim.min = NULL
	if (!is.numeric(ylim.min)) { # auto
		ylim.min = min(modkm$lower[!is.na(modkm$lower)]) - .2
		ylim.min = max(ylim.min, 0)
		params$ylim          = c(ylim.min, 1)
		params$pval.coord[2] = params$pval.coord[2] + ylim.min
	} else if (ylim.min > 0) {
		params$ylim          = c(ylim.min, 1)
		params$pval.coord[2] = params$pval.coord[2] + ylim.min
	}
	params$fit     = modkm
	params$data    = data

	if (length(varlvls) < length(params$legend.labs))
		params$legend.labs = varlvls

	ret.plot      = do.call(ggsurvplot, params)
	ret.plot$plot = apply.ggs(ret.plot$plot, mainggs)

	# hack the table
	if (params$risk.table) {
		tableggs.one = list(
			ylab             = list(params$legend.title),
			xlab             = list(paste0("Time (", outunit ,")")),
			scale_y_discrete = list(labels = rev(params$legend.labs))
		)
		tableggs.one   = update.list(tableggs.one, tableggs)
		ret.plot$table = apply.ggs(ret.plot$table, tableggs.one)
	}
	do.call(png, c(paste0(prefix, var, '.survival.png'), devpars))
	print(ret.plot)
	dev.off()

	ret.summary = matrix(ncol = 6, nrow = 0)
	colnames(ret.summary) = c('var', 'method', 'test', 'groups', 'df', 'pvalue')
	ret.summary = rbind(ret.summary, c(v, 'logrank', 'NA', paste(table(data[, var]), collapse = ','), 1, survp$pval))
	write.table(ret.summary, paste0(prefix, var, '.survival.txt'), row.names = F, col.names = T, sep = "\t", quote = F)
}

plot.survival(survdata, 'LCSIG')
if (gene != "") {
	plot.survival(survdata, paste('LCSIG', gene, sep = ':'))
}