library(survival)
library(survminer)
{{rimport}}('__init__.r', 'plot.r')

set.seed(8525)

# load parameters
infile      = {{i.infile | R}}
outfile     = {{o.outfile | R}}
outdir      = {{o.outdir | R}}
inunit      = {{args.inunit | R}}
outunit     = {{args.outunit | R}}
covfile     = {{args.covfile | R}}
nthread     = {{args.nthread | R}}
rnames      = {{args.inopts | lambda x: x.get('rnames', True) | R}}
combine     = {{args.combine | R}}
method      = {{args.method | R}}
devpars     = {{args.devpars | R}}
params.plot = {{args.params | R}}
pval        = {{args.pval | R}}
mainggs     = {{args.ggs | R}}
ngroups     = as.numeric({{args.ngroups | R}})
autogroup   = as.logical({{args.autogroup | R}})

if (pval == T) {
	pval = 'logrank'
}

tableggs      = mainggs$table
mainggs$table = NULL

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

if (is.list(params.plot) && length(params.plot) > 0) {
	params.default = list(
		pval             = '{method}\np = {pval}',
		risk.table       = T,
		#surv.median.line = 'hv',
		conf.int         = T
	)
	params.plot = update.list(params.default, params.plot)
}

# cut the data in to n parts using quantile
# in case the data has equal numbers, add very small random numbers to it.
ncut = function(dat, n, prefix = 'q') {
	indexes = 1:length(dat)
	x       = cut(
		indexes, 
		breaks         = quantile(indexes, probs = seq(0,1,1/n)),
		include.lowest = T,
		labels         = paste0(prefix, 1:n)
	)
	return (x[match(indexes, order(dat))])
}

# get the median value of all types
allTypeMedian = function(dat) {
	return(dat[median(1:length(dat))])
}

# apply survival analysis method (cox or km)
useCox = function(dat, varname) {
	cnames  = colnames(dat)
	var     = cnames[3]
	# all variables
	vars    = cnames[3:length(cnames)]
	# levels of the variable
	varlvls = levels(factor(dat[, var]))
	# apply the model
	fmula.str = paste('Surv(time, status) ~', paste(vars, collapse = ' + '))
	modcox    = coxph(as.formula(fmula.str), data = as.data.frame(data.matrix(dat)))

	# the results
	#   - summary: The result table of the main variable
	#   - cov:  The result of the covariates
	#   - plot: The survival plot
	ret = list(summary = matrix(ncol = 8, nrow = 0), cov = matrix(ncol = 7, nrow = 0), plot = NULL)
	colnames(ret$summary) = c('var', 'method', 'test', 'df', 'pvalue', 'hr', 'ci95_lower', 'ci95_upper')
	colnames(ret$cov) = c('mainvar', 'covar', 'hr', 'ci95_lower', 'ci95_upper', 'z', 'pvalue')

	s = summary(modcox)
	ret$summary = rbind(ret$summary, c(varname, 'logrank', s$sctest, s$conf.int[var, -2]))
	ret$summary = rbind(ret$summary, c(varname, 'waldtest', s$waldtest, s$conf.int[var, -2]))
	ret$summary = rbind(ret$summary, c(varname, 'likeratio', s$logtest, s$conf.int[var, -2]))

	for (v in vars) {
		ret$cov = rbind(ret$cov, c(varname, v, s$conf.int[v, -2], s$coefficients[v, c(4,5)]))
	}
	ret$summary = pretty.numbers(ret$summary, list(
		test..hr..ci95_lower..ci95_upper = '%.3f',
		pvalue = '%.2E'
	))
	ret$cov = pretty.numbers(ret$cov, list(
		hr..ci95_lower..ci95_upper..z = '%.3f',
		pvalue = '%.2E'
	))

	# the pvaldata used to replace the place holders in args.plot.params.pval
	if (!is.logical(pval)) {
		pvaldata = list(
			method = test.names[[pval]],
			pval   = sprintf('%.2E', s[[ test.coxes[[pval]] ]][3]),
			hr     = sprintf('%.3f', s$conf.int[var, 1]),
			ci95L  = sprintf('%.3f', s$conf.int[var, 3]),
			ci95U  = sprintf('%.3f', s$conf.int[var, 4])
		)
	}

	# plot the result out
	if (is.list(params.plot) && length(params.plot) > 0) {
		params.default = list(
			legend.title = paste0('Strata(', varname, ')'),
			legend.labs  = ifelse(length(varlvls)<=ngroups, varlvls, paste0('q', 1:ngroups)),
			xlab         = paste0("Time (", outunit ,")"),
			ylim.min     = 0.0,
			pval.coord   = c(0, 0.1)
		)
		params          = update.list(params.default, params.plot)
		
		params$legend.title = gsub('{var}', varname, params$legend.title, fixed = T)
		if (!is.logical(pval)) {
			params$pval = gsub('{method}', pvaldata$method, params$pval, fixed = T)
			params$pval = gsub('{pval}',   pvaldata$pval,   params$pval, fixed = T)
			params$pval = gsub('{hr}',     pvaldata$hr,     params$pval, fixed = T)
			params$pval = gsub('{ci95L}',  pvaldata$ci95L,  params$pval, fixed = T)
			params$pval = gsub('{ci95U}',  pvaldata$ci95U,  params$pval, fixed = T)
		} else {
			params$pval = FALSE
		}
		# use Kaplan Merier to plot
		km            = useKM(dat, varname, params)
		ret$plot      = km$plot
		ret$quantdata = km$quantdata

		ret$summary = cbind(ret$summary, groups = rep(km$groups, nrow(ret$summary)))
	}

	return(ret)
}


useKM = function(dat, varname, params = NULL) {
	# @description:
	#   Use Kaplan Merier to do survival analysis
	# @params:
	#   dat: The data frame, where:
	#     - 1st col is time
	#     - 2nd col is status (censoring)
	#     - 3rd col is the varialbe (var)
	#   varname: The variable name to display instead of var (3rd colname)
	#   params:  The params used for ggsurvplot. If provided, global plot$params will not be override

	# get the variable name
	cnames  = colnames(dat)
	var     = cnames[3]

	# get the levels of the values of the variable
	varlvls = levels(factor(dat[, var]))
	# how many levels are there?
	nvlvls = length(varlvls)

	fmula = as.formula(paste('Surv(time, status) ~', var))
	# the result: 
	#   - summary: The result table
	#   - plot:    The survival plot
	#   - quantdata: The data to determine the quantile percentage
	ret   = list(summary = matrix(ncol = 6, nrow = 0), plot = NULL, quantdata = NULL)

	# if there are way more levels than expected number of groups, do the cut
	if (nvlvls > ngroups) {
		if (ngroups == 2 && autogroup) {
			quantsteps = seq(0.15, 0.85, by = 0.005)
			survscores = sapply(quantsteps, function(x){
				newdata = dat
				newdata[, var] = as.numeric(dat[, var] >= quantile(dat[, var], x))
				tryCatch({
					coxph(fmula, data = newdata)$score
				}, error = function(x) 1)
			})
			if (is.list(params.plot)) {
				ret$quantdata = data.frame(step = quantsteps, survscore = survscores)
			}
			idxmax = sample(which(survscores == max(survscores)), 1)
			quant  = quantsteps[idxmax]
			dat[, var] = paste0(
				"q", 
				c(1, 2)[as.numeric(dat[, var] >= quantile(dat[, var], quant)) + 1]
			)
			varlvls = levels(factor(dat[, var]))
			rm(quantsteps, survscore)
		} else {
			dat[, var] = ncut(dat[, var], ngroups)
			varlvls = levels(factor(dat[, var]))
		}
	}
	# apply the model
	modkm     = surv_fit(fmula, data = dat)

	# called directly
	if (is.null(params)) {
		# get the pvalue
		survpval  = surv_pvalue(modkm)
		pval      = 'logrank'
		# format the pvalue
		if (!is.logical(pval)) {
			pvaldata  = list(
				method = test.names[[pval]],
				pval   = sprintf('%.2E', survpval$pval)
			)
		}

		colnames(ret$summary) = c('var', 'method', 'test', 'groups', 'df', 'pvalue')
		if (nvlvls == 1) {
			ret$summary = rbind(ret$summary, c(varname, 'logrank', 'NA', nrow(dat), 1, '1.00E+00'))
			return (ret)
		} else {
			ret$summary = rbind(ret$summary, c(varname, 'logrank', 'NA', paste(table(dat[, var]), collapse = ','), 1, survpval$pval))
			ret$summary = pretty.numbers(ret$summary, list(
				pvalue = '%.2E'
			))
		}
	}

	if (is.list(params.plot) && length(params.plot) > 0) {
		if (is.null(params)) {
			params.default = list(
				legend.title = varname,
				legend.labs  = varlvls,
				xlab         = paste0("Time (", outunit ,")"),
				ylim.min     = 0.0,
				pval.coord   = c(0, 0.1)
			)
			params          = update.list(params.default, params.plot)
			
			if (!is.logical(pval)) {
				params$pval = gsub('{method}', pvaldata$method, params$pval, fixed = T)
				params$pval = gsub('{pval}',   pvaldata$pval,   params$pval, fixed = T)
			} else {
				params$pval = FALSE
			}
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
		params$data    = dat

		varlvls = levels(factor(dat[, var]))
		if (length(varlvls) < length(params$legend.labs))
			params$legend.labs = varlvls

		ret$groups    = paste(table(dat[, var]), collapse = ',')
		ret$plot      = do.call(ggsurvplot, params)
		ret$plot$plot = apply.ggs(ret$plot$plot, mainggs)

		# hack the table
		if (params$risk.table) {
			tableggs.one = list(
				ylab             = list(varname),
				xlab             = list(paste0("Time (", outunit ,")")),
				scale_y_discrete = list(labels = rev(params$legend.labs))
			)
			tableggs.one   = update.list(tableggs.one, tableggs)
			ret$plot$table = apply.ggs(ret$plot$table, tableggs.one)
		}
	}
	return(ret)
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
survivalOne = function(data, varname, cvdata = NULL) {
	if (!is.null(cvdata)) {
		data = cbind.fill(data, cvdata)
	}
	# only keep samples with all data available
	data      = data[complete.cases(data), , drop = F]
	# construct formula like 'Surv(time, status) ~ age + sex' using column names
	cnames    = make.names(colnames(data))
	colnames(data) = cnames
	
	var     = cnames[3]
	vars    = cnames[3:length(cnames)]
	varlvls = levels(factor(data[, var]))
	nvlvls  = length(varlvls)
	m       = method
	if (m == 'auto') {
		m = if (length(vars) == 1) 'km' else 'cox'
	}

	if (m == 'km') {
		ret = useKM(data, varname)
	} else {
		ret = useCox(data, varname)
	}

	return(list(ret))
}

data  = read.table.nodup(infile, sep = "\t", header = T, row.names = if(rnames) 1 else NULL, check.names = F)
vdata = NULL
if (!is.null(covfile) && covfile != '') {
	if (!rnames) 
		stop('Rownames are required for covariant analysis.')
	vdata = read.table.nodup(covfile, sep = "\t", header = T, row.names = 1, check.names = F)
}

fct = 1
{% case (args.inunit, args.outunit) %}
	{% when ('days', 'weeks') %}
	fct = 1 / 7
	{% when ('days' ,'months') %}
	fct = 1 / 30
	{% when ('days' ,'years') %}
	fct = 1 / 365
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
		rets = c(rets, survivalOne(data[, c(1, 2, i+2), drop = F], var, vdata))
	}
} else {
	library(doParallel)
	cl <- makeCluster(nthread)
	registerDoParallel(cl)

	rets = foreach (i = 1:length(vars), .verbose = T, .combine = c, .packages = c('survival', 'survminer')) %dopar% {
		survivalOne(data[, c(1, 2, i+2), drop = F], vars[i], vdata)
	}
	stopCluster(cl)
}


allsums = NULL
allcovs = NULL
if (is.list(combine) && length(combine) > 0) {
	survplots   = NULL
	for (ret in rets) {
		allsums = rbind(allsums, ret$summary)
		if ('cov' %in% names(ret))
			allcovs = rbind(allcovs, ret$cov)
		if (!is.null(ret$plot))
			survplots  = c(survplots, list(ret$plot))
		if (!is.null(ret$quantdata)) {
			plotfile = file.path(outdir, paste0(ret$summary[1, 'var'], '.quant.png'))
			plot.scatter(ret$quantdata, plotfile, ggs = list(geom_line = list()), devpars = devpars)
		}
	}
	if (!is.null(survplots)) {
		if (!'ncol' %in% names(combine)) combine$ncol = 1
		if (!'nrow' %in% names(combine)) combine$nrow = 1
		maxdim         = max(combine$ncol, combine$nrow)
		devpars$height = maxdim * devpars$height
		devpars$width  = maxdim * devpars$width
		plotfile = file.path(outdir, '{{i.infile | fn2}}.survival.png')
		do.call(png, c(plotfile, devpars))
		do.call(arrange_ggsurvplots, c(list(x = survplots, print = T), combine))
		dev.off()
	}
} else {
	for (ret in rets) {
		if ('cov' %in% names(ret))
			allcovs = rbind(allcovs, ret$cov)
		allsums  = rbind(allsums, ret$summary)
		if (!is.null(ret$plot)) {
			plotfile = file.path(outdir, paste0(ret$summary[1, 'var'], '.survival.png'))
			do.call(png, c(plotfile, devpars))
			print (ret$plot)
			dev.off()
		}
		if (!is.null(ret$quantdata)) {
			plotfile = file.path(outdir, paste0(ret$summary[1, 'var'], '.quant.png'))
			plot.scatter(ret$quantdata, plotfile, ggs = list(geom_line = list()), devpars = devpars)
		}
	}
}
write.table(allsums, outfile, sep="\t", quote=F, col.names = T, row.names = F)

if (!is.null(allcovs))
	write.table(allcovs, '{{o.outfile | prefix | prefix}}.covariants.txt', sep="\t", quote=F, col.names = T, row.names = F)
