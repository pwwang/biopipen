library(survival)
library(survminer)
{{rimport}}('__init__.r', 'plot.r')

set.seed(8525)

# load parameters
infile      = {{i.infile | R}}
incovfile   = {{i.covfile | R}}
outfile     = {{o.outfile | R}}
outdir      = {{o.outdir | R}}
inunit      = {{args.inunit | R}}
outunit     = {{args.outunit | R}}
inopts      = {{args.inopts | R}}
argscovfile = {{args.covfile | R}}
combine     = {{args.combine | R}}
method      = {{args.method | R}}
devpars     = {{args.devpars | R}}
params      = {{args.params | R}}
ggs         = {{args.ggs | R}}
covfile = ifelse(is.true(incovfile), incovfile, argscovfile)

data  = read.table.inopts(infile, inopts)

fct = 1
{% case (args.inunit, args.outunit) %}
	{% when ('days', 'weeks') %}
	fct = 1 / 7
	{% when ('days' ,'months') %}
	fct = 1 / 30
	{% when ('days' ,'years') %}
	fct = 1 / 365.25
	{% when ('weeks', 'days') %}
	fct = 7
	{% when ('weeks', 'months') %}
	fct = 7*12/365.25
	{% when ('weeks', 'years') %}
	fct = 7/365.25
	{% when ('years', 'days') %}
	fct = 365.25
	{% when ('years', 'weeks') %}
	fct = 365.25/7
	{% when ('years', 'months') %}
	fct = 12
{% endcase %}

if (fct != 1) {
	data[,1] = data[,1]*fct
}

cnames = make.names(colnames(data))
vars   = cnames[3:length(cnames)]
cnames[1] = 'Time'
cnames[2] = 'Status'

colnames(data) = cnames
covdata = NULL
if (is.true(covfile)) {
	if (method == 'km') {
		stop('Kaplan Merier is not able to do multi-variate analysis.')
	}
	covdata = read.table.inopts(covfile, list(rnames = TRUE, cnames = TRUE))
	rnames  = intersect(rownames(data), rownames(covdata))
	data    = data[rnames,,drop = FALSE]
	covdata = covdata[rnames,,drop = FALSE]
}

results = lapply(vars, function(v) {
	plotfile = file.path(outdir, sprintf('%s.surv.png', v))
	if (method == 'km') {
		r = plot.survival.km(data[, c(cnames[1:2], v), drop = FALSE],
			plotfile = plotfile, params = params, ggs = ggs, devpars = devpars)
	} else {
		if (!is.null(covdata)) {
			coxdata = cbind(data[, c(cnames[1:2], v), drop = FALSE], covdata)
		} else {
			coxdata = data[, c(cnames[1:2], v), drop = FALSE]
		}
		coxdata = coxdata[complete.cases(coxdata),,drop = FALSE]
		r = plot.survival.cox(coxdata,
			plotfile = plotfile, params = params, ggs = ggs, devpars = devpars)
	}
	as.list(r)
})

retable = NULL
# combine the tables
for (i in 1:length(results)) {
	if (method == 'cox') {
		results[[i]]$table = cbind(mainvar = vars[i], results[[i]]$table)
	}
	retable = ifelse(is.null(retable), results[[i]]$table, rbind(retable, results[[i]]$table))
}
write.table(pretty.numbers2(retable, HR..CI95_1..CI95_2..z = '%.3f', pvalue..modpval = '%.3E'), outfile, sep = "\t", row.names = FALSE, quote = FALSE)

if (is.true(combine) && length(results) > 1) {
	if(is.true(combine$nrow)) {
		devpars$height = devpars$height*combine$nrow
	}
	if(is.true(combine$ncol)) {
		devpars$width = devpars$width*combine$ncol
	}
	library(ggpubr)
	if (is.true(params$risk.table)) {

	}
	p = do.call(arrange_ggsurvplots,
		c(list(x=lapply(1:length(results), function(i) results[[i]]$plot), print=TRUE), combine))
	save.plot(p, file.path(outdir, 'survival.combined.png'), devpars)
}