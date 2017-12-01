library(survival)
library(survminer)

survivalOne = function(dat, var) {
	fmula = paste0('Surv({{args.outunit}}, status) ~ ', var)
	fmula = as.formula(fmula)
	sfit  = do.call(survfit, list(formula = fmula, data = dat))
	diff  = pairwise_survdiff(fmula, data=dat)
	#pvals = diff$p.value
	rns   = rownames(diff$p.value)
	cns   = colnames(diff$p.value)
	lenrs = length(rns)
	lencs = length(cns)
	pvals = NULL
	for (i in 1:lenrs) {
		for (j in 1:lencs) {
			if (i < j) next
			pvals = rbind(pvals, c(var = var, group1 = rns[i], group2 = cns[j], pval = formatC(diff$p.value[rns[i], cns[j]], format='e', digits = 2)))
		}
	}
	psurv = do.call(ggsurvplot, c(list(fit = sfit, pval = {{args.pval | R}}, legend.title = var, legend.labs = levels(dat[,var]), xlab = "Time ({{args.outunit}})"), {{args.plotParams | Rlist}}))
	return(list(ret = list(pvals = pvals, psurv = psurv, var = var)))
}

rnames = {% if args.rnames %}1{% else %}NULL{% endif %}
df = read.table({{in.infile | quote}}, sep = "\t", header = T, row.names = rnames, check.names = F)

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

cnames = colnames(df)
cnames[1] = {{args.outunit | quote}}
cnames[2] = 'status'
colnames(df) = cnames
if (fct != 1) {
	df[,1] = df[,1]*fct
}

vars = cnames[3:length(cnames)]
if (length(vars) == 1 || {{args.nthread}} == 1) {
	rets = NULL
	for (var in vars)  {
		rets  = c(rets, survivalOne(df, var))
	}
} else {
	library(doParallel)
	cl <- makeCluster({{args.nthread}})  
    registerDoParallel(cl)

	rets = foreach (i = 1:length(vars), .verbose = T, .combine = c, .packages = c('survival', 'survminer')) %dopar% {
		survivalOne(df, vars[i])
	}
	stopCluster(cl) 
}

devpars    = {{args.devpars | Rlist}}
gridParams = {{args.gridParams | Rlist}}
if (!'ncol' %in% names(gridParams)) gridParams$ncol = 1
if (!'nrow' %in% names(gridParams)) gridParams$nrow = 1
maxdim         = max(gridParams$ncol, gridParams$nrow)
devpars$height = maxdim * devpars$height
devpars$width  = maxdim * devpars$width

outpval = file.path({{out.outdir | quote}}, 'survival.pval.txt') 
pvals   = NULL
{% if args.combine %}
outplot = file.path({{out.outdir | quote}}, 'survival.png')
psurvs  = NULL
for (ret in rets) {
	pvals = rbind(pvals, ret$pvals)
	psurvs = c(psurvs, list(ret$psurv))
}
do.call (png, c(outplot, devpars))
do.call (arrange_ggsurvplots, c(list(x = psurvs, print = T), gridParams))
dev.off()
{% else %}
for (ret in rets) {
	pvals   = rbind(pvals, ret$pvals)
	outpv   = file.path({{out.outdir | quote}}, paste0(ret$var, '.pval.txt'))
	outplot = file.path({{out.outdir | quote}}, paste0(ret$var, '.png'))
	write.table(ret$pvals, outpv, sep="\t", quote=F, col.names = T, row.names = F)
	do.call (png, c(outplot, devpars))
	print (ret$psurv)
	dev.off()
}
{% endif %}
write.table(pvals, outpval, sep="\t", quote=F, col.names = T, row.names = F)
