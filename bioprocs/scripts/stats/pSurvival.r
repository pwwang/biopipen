library(survival)
library(survminer)
{{rimport}}('__init__.r', 'plot.r')

# load parameters
infile  = {{in.infile | R}}
outfile = {{out.outfile | R}}
outdir  = {{out.outdir | R}}
inunit  = {{args.inunit | R}}
outunit = {{args.outunit | R}}
covfile = {{args.covfile | R}}
nthread = {{args.nthread | R}}
rnames  = {{args.inopts | lambda x: x.get('rnames', True) | R}}
combine = {{args.combine | R}}
devpars = {{args.devpars | R}}
plot    = {{args.plot | R}}
pval    = {{args.pval | R}}
ggs     = {{args.ggs | R}}

if (pval == T) {
	pval = 'logrank'
}

test.names = list(
	logrank   = 'Logrank test',
	waldtest  = 'Walk test',
	likeratio = 'Likelihood ratio test'
)
test.coxes = list(
	logrank   = 'sctest',
	waldtest  = 'waldtest',
	likeratio = 'logtest'
)

# partition the variable into n-tiles
partVar = function(data, n, prefix = 'q') {
	tiles = quantile(data, probs = seq(0, 1, 1/n))
	return(unlist(lapply(
		data, 
		function(x) {
			for (i in n:1) {
				if (x >= tiles[i]) return(paste0(prefix, i))
			}
		}
	)))
}

# do survival analysis on one variable and all covariants
# the input data should be a data frame with header:
# 
# 	time	status	var
# sample1	100	1	1.2
# sample2	89	0	3.3
# 
# the cvdata (covariant data frame) should be:
# 	sex	age
# sample1	1	39
# sample2	0	87
survivalOne = function(data, cvdata = NULL) {
	if (!is.null(cvdata)) {
		data = cbind.fill(data, cvdata)
	}
	# only keep samples with all data available
	data = data[complete.cases(data), , drop = F]
	# construct formula like 'Surv(time, status) ~ age + sex' using column names
	vars      = colnames(data)
	vars      = vars[3:length(vars)]
	var       = vars[1]
	fmula.str = paste('Surv(time, status)', '~', paste(vars, collapse = ' + '))
	fmula     = as.formula(fmula.str)
	coxmodel  = coxph(fmula, data = data)
	coxret    = summary(coxmodel)

	# plot: The ggsurvplot object
	# test: The test result from different tests
	ret = list(plot = NULL, test = matrix(ncol = 5, nrow = 0))
	colnames(ret$test) = c('var', 'method', 'test', 'df', 'pvalue')

	ret$test = rbind(ret$test, c(var, 'logrank', coxret$sctest))
	ret$test = rbind(ret$test, c(var, 'waldtest', coxret$waldtest))
	ret$test = rbind(ret$test, c(var, 'likeratio', coxret$logtest))

	if ((class(plot) == 'list' && length(plot) == 0) || (class(plot) == 'logical' && !plot))
		return(list(ret))
	
	# construct new data to plot according to args.plot.var
	varlvls = levels(factor(data[, var]))
	ncurves = plot$ncurves
	if (length(varlvls) <= ncurves) {
		ncurves = length(varlvls)
		vardata = data[, var]
	} else {
		vardata = partVar(data[, var, drop = T], ncurves)
	}
	newdata = data.frame(var = levels(factor(vardata)))
	colnames(newdata) = var
	for (v in vars) {
		if (v == var) next
		tmpdata = matrix(rep(mean(data[, v], na.rm = TRUE), ncurves), ncol = 1)
		colnames(tmpdata) = v
		newdata = cbind(newdata, tmpdata)
	}

	varcount = table(vardata)
	
	fit         = survfit(coxmodel, newdata = newdata)
	params      = list(
		risk.table   = T,
		pval         = '{method}\np = {pval}',
		legend.title = paste0('Strata(', var, ')'),
		legend.labs  = paste0(newdata[, var], '(', varcount, ')'),
		xlab         = paste0("Time (", outunit ,")")
	)
	params$fit  = fit
	params$data = newdata
	params      = update.list(params, plot$params)
	if (pval != F) {
		params$pval = gsub('{method}', test.names[[pval]], params$pval, fixed = T)
		params$pval = gsub('{pval}', sprintf('%.2E', coxret[[ test.coxes[[pval]] ]][3]), params$pval, fixed = T)
	}
	ret$plot      = do.call(ggsurvplot, params)
	mainggs       = ggs
	mainggs$table = NULL
	ret$plot$plot = apply.ggs(ret$plot$plot, mainggs)
	if (params$risk.table) {
		kmdata               = data
		kmdata[, var]        = vardata
		kmfmula.str          = paste('Surv(time, status)', '~', var)
		kmfmula              = as.formula(kmfmula.str)
		kmfit                = surv_fit(kmfmula, data = kmdata)
		kmparams             = params
		kmparams$fit         = kmfit
		kmparams$data        = kmdata
		kmparams$risk.table  = T
		kmparams$legend.labs = newdata[, var]
		kmsurv               = do.call(ggsurvplot, kmparams)
		ret$plot$table       = apply.ggs(kmsurv$table, ggs$table)
	}

	return(list(ret))
}

data  = read.table.nodup(infile, sep = "\t", header = T, row.names = if(rnames) 1 else NULL, check.names = F)
vdata = NULL
if (!is.null(covfile)) {
	if (!rnames) 
		stop('Rownames are required for covariant analysis.')
	vdata = read.table.nodup(covfile, sep = "\t", header = T, row.names = 1, check.names = F)
}

fct = 1
{% if args.inunit, args.outunit | lambda x,y: x == 'days' and y == 'weeks' %}
fct = 1 / 7
{% elif args.inunit, args.outunit | lambda x,y: x == 'days' and y == 'months' %}
fct = 1 / 30
{% elif args.inunit, args.outunit | lambda x,y: x == 'days' and y == 'years' %}
fct = 1 / 365
{% elif args.inunit, args.outunit | lambda x,y: x == 'weeks' and y == 'days' %}
fct = 7
{% elif args.inunit, args.outunit | lambda x,y: x == 'weeks' and y == 'months' %}
fct = 7*12/365
{% elif args.inunit, args.outunit | lambda x,y: x == 'weeks' and y == 'years' %}
fct = 7/365
{% elif args.inunit, args.outunit | lambda x,y: x == 'years' and y == 'days' %}
fct = 365
{% elif args.inunit, args.outunit | lambda x,y: x == 'years' and y == 'weeks' %}
fct = 365/7
{% elif args.inunit, args.outunit | lambda x,y: x == 'years' and y == 'months' %}
fct = 12
{% endif %}

cnames = colnames(data)
vars   = cnames[3:length(cnames)]

cnames[1] = 'time'
cnames[2] = 'status'

colnames(data) = cnames

if (fct != 1) {
	data[,1] = data[,1]*fct
}

# if there is only one variable to analyze
if (length(vars) == 1 || nthread == 1) {
	rets = NULL
	for (i in 1:length(vars))  {
		var  = vars[i]

		rets = c(rets, survivalOne(data[, c(1, 2, i+2), drop = F], vdata))
	}
} else {
	library(doParallel)
	cl <- makeCluster(nthread)
	registerDoParallel(cl)

	rets = foreach (i = 1:length(vars), .verbose = T, .combine = c, .packages = c('survival', 'survminer')) %dopar% {
		survivalOne(data[, c(1, 2, i+2), drop = F], vdata)
	}
	stopCluster(cl)
}


allouts = NULL
if (combine) {
	psurvs   = NULL
	for (ret in rets) {
		allouts = rbind(allouts, ret$test)
		if (!is.null(ret$plot))
			psurvs  = c(psurvs, list(ret$plot))
	}
	if (!is.null(psurvs)) {
		if (!'ncol' %in% names(plot$arrange)) plot$arrange$ncol = 1
		if (!'nrow' %in% names(plot$arrange)) plot$arrange$nrow = 1
		maxdim         = max(plot$arrange$ncol, plot$arrange$nrow)
		devpars$height = maxdim * devpars$height
		devpars$width  = maxdim * devpars$width
		plotfile = file.path(outdir, '{{in.infile | fn2}}.survival.png')
		do.call(png, c(plotfile, devpars))
		do.call(arrange_ggsurvplots, c(list(x = psurvs, print = T), plot$arrange))
		dev.off()
	}
} else {
	for (ret in rets) {
		allouts  = rbind(allouts, ret$test)
		if (!is.null(ret$plot)) {
			plotfile = file.path(outdir, paste0(ret$test[1, 'var'], '.survival.png'))
			do.call(png, c(plotfile, devpars))
			print (ret$plot)
			dev.off()
		}
	}
}
write.table(pretty.numbers(allouts, formats = list(pvalue = '%.2E', test = '%.3f')), outfile, sep="\t", quote=F, col.names = T, row.names = F)
