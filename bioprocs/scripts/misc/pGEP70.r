{{rimport}}('__init__.r', 'plot.r')

{% python from os import path %}
library(methods)
library(survival)
library(survminer)

inunit   = {{args.inunit | R}}
outunit  = {{args.outunit | R}}
exprfile = {{i.exprfile | R}}
survfile = {{i.survfile | R}}
gene     = {{i.gene | quote}}
gep70    = {{args.gep70 | R}}
devpars  = {{args.devpars | R}}
params   = {{args.params | R}}
name     = {{args.name | R}}
prefix   = {{fn2(i.survfile), i.gene, o.outdir, args.name | lambda v1, v2, v3, v4, path = path: path.join(v3, '{}.{}-{}.'.format(v1, v4, v2) if v2 else '{}.{}.'.format(v1, v4)) | quote }}

mainggs       = {{args.ggs | R}}
tableggs      = mainggs$table
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
if (fct != 1) {
	survdata[, 1] = survdata[, 1] * fct
}

# read all expression data
# rows: genes
# cols: samples
exprdata = read.table.nodup(exprfile, header = T, row.names = 1, sep = "\t", check.names = F)
allgenes = rownames(exprdata)

# read the 70 genes
gene70 = read.table(gep70, row.names = 1, header = T, sep = "\t", check.names = F)

gene51   = intersect(allgenes, rownames(gene70[which(gene70[,1] == 'up'),,drop=F]))
gene19   = intersect(allgenes, rownames(gene70[which(gene70[,1] == 'down'),,drop=F]))	
exp51    = as.matrix(colMeans(exprdata[gene51, , drop = F]))
exp19    = as.matrix(colMeans(exprdata[gene19, , drop = F]))
gep70val = exp51 - exp19
colnames(gep70val) = name

samples  = intersect(rownames(gep70val), rownames(survdata))
gep70val = gep70val[samples, , drop=F]
if (gene != "") {
	gexpr = exprdata[gene, samples, drop = F]
	survdata = cbind(survdata[samples, ], gep70val, t(gexpr))
} else {
	survdata = cbind(survdata[samples, ], gep70val)
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

survdata[, name] = c('low', 'high')[as.numeric(survdata[, name] > quantile(
	survdata[, name], 
	autogroup(survdata, name)
)) + 1]

if (gene != "") {
	survdata[, gene] = c('low', 'high')[as.numeric(survdata[, gene] > quantile(
		survdata[, gene], 
		autogroup(survdata, gene)
	)) + 1]
	Group = apply(survdata, 1, function(row) paste(row[3], row[4], sep = ':'))
	Group = matrix(Group, ncol = 1)
	colnames(Group) = paste(name, gene, sep = ':')
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

plot.survival(survdata, name)
if (gene != "") {
	plot.survival(survdata, paste(name, gene, sep = ':'))
}

