library(survival)
library(survminer)
{{rimport}}('__init__.r', 'plot.r')

set.seed(8525)

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
method  = {{args.method | R}}
devpars = {{args.devpars | R}}
plot    = {{args.plot | R}}
pval    = {{args.pval | R}}
mainggs = {{args.ggs | R}}
ngroups = {{args.ngroups | R}}

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

if (is.list(plot)) {
	params = list(
		pval             = '{method}\np = {pval}',
		risk.table       = T,
		surv.median.line = 'hv',
		conf.int         = T
	)
	plot$params = update.list(params, plot$params)
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

# convert group names to numeric
to.numeric = function(dat) {
	return(as.numeric(gsub('[^0-9.+-]', '', dat)))
}

# apply survival analysis method (cox or km)
useCox = function(dat, varname) {
	cnames  = colnames(dat)
	var     = cnames[3]
	vars    = cnames[3:length(cnames)]
	varlvls = levels(factor(dat[, var]))
	if (!is.numeric(varlvls)) {
		dat[, var] = to.numeric(dat[, var])
	}
	fmula.str = paste('Surv(time, status) ~', paste(vars, collapse = ' + '))
	modcox    = coxph(as.formula(fmula.str), data = dat)

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
	if (!is.logical(pval)) {
		pvaldata = list(
			method = test.names[[pval]],
			pval   = sprintf('%.2E', s[[ test.coxes[[pval]] ]][3]),
			hr     = sprintf('%.3f', s$conf.int[var, 1]),
			ci95L  = sprintf('%.3f', s$conf.int[var, 3]),
			ci95U  = sprintf('%.3f', s$conf.int[var, 4])
		)
	}

	if (is.list(plot)) {
		nvlvls  = length(varlvls)
		kmdata = dat
		if (nvlvls > ngroups) {
			kmdata[, var] = ncut(kmdata[, var], ngroups)
		}
		kmfmula.str = paste('Surv(time, status) ~', var)
		kmfmula     = as.formula(kmfmula.str)
		# survfit leads to an error
		# see: https://github.com/kassambara/survminer/issues/283
		kmfit       = surv_fit(kmfmula, data = kmdata)

		if (nvlvls <= ngroups) {
			newdata = data.frame(.var = levels(factor(dat[, var])))
			colnames(newdata) = var
			for (v in vars) {
				if (var == v) next
				tmp = data.frame(.v = allTypeMedian(dat[, v]))
				colnames(tmp) = v
				newdata = cbind(newdata, tmp)
			}
		} else {
			newvars = ncut(dat[, var], ngroups)
			varlvls = levels(factor(newvars))
			newdata = data.frame(.var = unlist(lapply(varlvls, function(x) mean(dat[which(newvars == x), var]))))
			colnames(newdata) = var
			for (v in vars) {
				if (var == v) next
				tmp = data.frame(.v = allTypeMedian(dat[, v]))
				colnames(tmp) = v
				newdata = cbind(newdata, tmp)
			}
		}
		params.default = list(
			legend.title = paste0('Strata(', varname, ')'),
			legend.labs  = varlvls,
			xlab         = paste0("Time (", outunit ,")"),
			ylim.min     = 0.0,
			pval.coord   = c(0, 0.1)
		)
		params          = update.list(params.default, plot$params)
		ylim.min        = params$ylim.min
		params$ylim.min = NULL
		if (!is.numeric(ylim.min)) { # auto
			ylim.min = min(kmfit$lower) - .2
			ylim.min = max(ylim.min, 0)
			params$ylim          = c(ylim.min, 1)
			params$pval.coord[2] = params$pval.coord[2] + ylim.min
		} else if (ylim.min > 0) {
			params$ylim          = c(ylim.min, 1)
			params$pval.coord[2] = params$pval.coord[2] + ylim.min
		}
		params$fit     = survfit(modcox, newdata = newdata)
		params$data    = newdata
		if (!is.logical(pval)) {
			params$pval = gsub('{method}', pvaldata$method, params$pval, fixed = T)
			params$pval = gsub('{pval}',   pvaldata$pval,   params$pval, fixed = T)
			params$pval = gsub('{hr}',     pvaldata$hr,     params$pval, fixed = T)
			params$pval = gsub('{ci95L}',  pvaldata$ci95L,  params$pval, fixed = T)
			params$pval = gsub('{ci95U}',  pvaldata$ci95U,  params$pval, fixed = T)
		} else {
			params$pval = FALSE
		}
		ret$plot       = do.call(ggsurvplot, params)
		ret$plot$plot  = apply.ggs(ret$plot$plot, mainggs)

		# hack the table
		if (params$risk.table) {
			tableggs.one = list(
				xlab             = list(paste0("Time (", outunit ,")")),
				ylab             = list(paste0('Strata(', varname, ')')),
				scale_y_discrete = list(labels = rev(varlvls))
			)
			tableggs.one   = update.list(tableggs.one, tableggs)
			ret$plot$table = apply.ggs(ggsurvplot(kmfit, data = kmdata, risk.table = TRUE)$table, tableggs.one)
		}
	}

	return(ret)
}

useKM = function(dat, varname) {
	cnames  = colnames(dat)
	var     = cnames[3]
	varlvls = levels(factor(dat[, var]))
	nvlvls  = length(varlvls)
	if (nvlvls > ngroups) {
		dat[, var] = ncut(dat[, var], ngroups)
		varlvls = levels(factor(dat[, var]))
	}
	fmula.str = paste('Surv(time, status) ~', var)
	modkm     = surv_fit(as.formula(fmula.str), data = dat)

	survpval  = surv_pvalue(modkm)
	pval      = 'logrank'
	if (!is.logical(pval)) {
		pvaldata  = list(
			method = test.names[[pval]],
			pval   = sprintf('%.2E', survpval$pval)
		)
	}

	ret = list(summary = matrix(ncol = 5, nrow = 0), plot = NULL)
	colnames(ret$summary) = c('var', 'method', 'test', 'df', 'pvalue')
	if (nvlvls == 1) {
		ret$summary = rbind(ret$summary, c(varname, 'logrank', 'NA', 1, '1.00E+00'))
		return (ret)
	} else {
		ret$summary = rbind(ret$summary, c(varname, 'logrank', 'NA', 1, survpval$pval))
		ret$summary = pretty.numbers(ret$summary, list(
			pvalue = '%.2E'
		))
	}

	if (is.list(plot)) {
		params.default = list(
			legend.title = varname,
			legend.labs  = varlvls,
			xlab         = paste0("Time (", outunit ,")"),
			ylim.min     = 0.0,
			pval.coord   = c(0, 0.1)
		)
		params          = update.list(params.default, plot$params)
		ylim.min        = params$ylim.min
		params$ylim.min = NULL
		if (!is.numeric(ylim.min)) { # auto
			ylim.min = min(modkm$lower) - .2
			ylim.min = max(ylim.min, 0)
			params$ylim          = c(ylim.min, 1)
			params$pval.coord[2] = params$pval.coord[2] + ylim.min
		} else if (ylim.min > 0) {
			params$ylim          = c(ylim.min, 1)
			params$pval.coord[2] = params$pval.coord[2] + ylim.min
		}
		params$fit     = modkm
		params$data    = dat
		if (!is.logical(pval)) {
			params$pval = gsub('{method}', pvaldata$method, params$pval, fixed = T)
			params$pval = gsub('{pval}',   pvaldata$pval,   params$pval, fixed = T)
		} else {
			params$pval = FALSE
		}
		ret$plot       = do.call(ggsurvplot, params)
		ret$plot$plot  = apply.ggs(ret$plot$plot, mainggs)

		# hack the table
		if (params$risk.table) {
			tableggs.one = list(
				ylab             = list(varname),
				xlab             = list(paste0("Time (", outunit ,")")),
				scale_y_discrete = list(labels = rev(varlvls))
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
if (combine) {
	survplots   = NULL
	for (ret in rets) {
		allsums = rbind(allsums, ret$summary)
		if ('cov' %in% names(ret))
			allcovs = rbind(allcovs, ret$cov)
		if (!is.null(ret$plot))
			survplots  = c(survplots, list(ret$plot))
	}
	if (!is.null(survplots)) {
		if (!'ncol' %in% names(plot$arrange)) plot$arrange$ncol = 1
		if (!'nrow' %in% names(plot$arrange)) plot$arrange$nrow = 1
		maxdim         = max(plot$arrange$ncol, plot$arrange$nrow)
		devpars$height = maxdim * devpars$height
		devpars$width  = maxdim * devpars$width
		plotfile = file.path(outdir, '{{in.infile | fn2}}.survival.png')
		do.call(png, c(plotfile, devpars))
		do.call(arrange_ggsurvplots, c(list(x = survplots, print = T), plot$arrange))
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
	}
}
write.table(allsums, outfile, sep="\t", quote=F, col.names = T, row.names = F)

if (!is.null(allcovs))
	write.table(allcovs, '{{out.outfile | prefix | prefix}}.covariants.txt', sep="\t", quote=F, col.names = T, row.names = F)
