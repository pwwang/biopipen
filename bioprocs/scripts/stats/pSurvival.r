library(survival)
library(survminer)

survivalOne = function(dat, var, vname) {
	# remove NAs
	dat  = dat[!is.na(dat[,var]), ,drop=F]
	lvls = levels(factor(dat[,var]))
	if (length(lvls) == 1) {
		cat(paste0("pyppl.log.warnning: Cannot do survival analysis for '",var,"' with one level.\n"), file=stderr())
		return(NULL)
	}
	freqs = as.data.frame(table(dat[,var]))
	fmula = paste0('Surv({{args.outunit}}, status) ~ ', var)
	fmula = as.formula(fmula)
	sfit  = do.call(survfit, list(formula = fmula, data = dat))
	diff  = pairwise_survdiff(fmula, data=dat)
	
	rns   = rownames(diff$p.value)
	if (length(rns) == 0) {
		pvals   = data.frame(var = vname, group1 = "", group2 = "", pval = 1, n1 = "", n2 = "")
		globalp = 1
	} else {
		cns   = colnames(diff$p.value)
		lenrs = length(rns)
		lencs = length(cns)
		pvals = NULL
		for (i in 1:lenrs) {
			for (j in 1:lencs) {
				if (i < j) next
				pvals = rbind(pvals, c(
					var = vname, 
					group1 = rns[i], 
					group2 = cns[j], 
					pval = formatC(diff$p.value[rns[i], cns[j]], format='e', digits = 2),
					n1 = freqs[freqs[,1] == rns[i], "Freq"],
					n2 = freqs[freqs[,1] == rns[j], "Freq"]
				))
			}
		}

		globaldiff = survdiff(fmula, data = dat)
		globalp    = 1 - pchisq(globaldiff$chisq, length(globaldiff$n) - 1)
	}
	labs = unlist(lapply(lvls, function(x) paste0(x, ' (n=', freqs[freqs[,1] == x, "Freq"], ')')))
	psurv = do.call(ggsurvplot, c(list(fit = sfit, pval = {{args.pval | R}}, legend.title = vname, legend.labs = labs, xlab = "Time ({{args.outunit}})"), {{args.plotParams | Rlist}}))
	return(list(ret = list(pvals = pvals, psurv = psurv, var = vname, globalp = globalp)))
}

rnames = {% if args.rnames %}1{% else %}NULL{% endif %}
df = read.table({{in.infile | quote}}, sep = "\t", header = T, row.names = rnames, check.names = F)

if (nrow(df) == 0 && {{args.noerror | R}}) {
	cat ('pyppl.log.warning: No survival data found.\n', file=stderr())
	writeLines("var\tglobalp\n", file.path({{out.outdir | quote}}, 'global.pval.txt'))
} else {

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
	vnames = cnames[3:length(cnames)]
	cnames    = make.names(cnames)
	cnames[1] = {{args.outunit | quote}}
	cnames[2] = 'status'
	colnames(df) = cnames
	if (fct != 1) {
		df[,1] = df[,1]*fct
	}
	
	survdat = df[, c(1,2), drop=F]
	vars = cnames[3:length(cnames)]
	if (length(vars) == 1 || {{args.nthread}} == 1) {
		rets = NULL
		for (i in 1:length(vars))  {
			rets  = c(rets, survivalOne(cbind(survdat, df[, vars[i], drop=F]), vars[i], vnames[i]))
		}
	} else {
		library(doParallel)
		cl <- makeCluster({{args.nthread}})  
		registerDoParallel(cl)

		rets = foreach (i = 1:length(vars), .verbose = T, .combine = c, .export = "survivalOne", .packages = c('survival', 'survminer')) %dopar% {
			survivalOne(cbind(survdat, df[, vars[i], drop=F]), vars[i], vnames[i])
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

	glbpval = file.path({{out.outdir | quote}}, 'global.pval.txt') 
	outpval = file.path({{out.outdir | quote}}, 'survival.pval.txt') 
	pvals   = NULL
	glbpv   = matrix(, ncol=2, nrow = 0)
	{% if args.combine %}
	outplot = file.path({{out.outdir | quote}}, 'survival.png')
	psurvs  = NULL
	for (ret in rets) {
		pvals = rbind(pvals, ret$pvals)
		glbpv = rbind(glbpv, c(ret$var, ret$globalp))
		psurvs = c(psurvs, list(ret$psurv))
	}
	do.call (png, c(outplot, devpars))
	do.call (arrange_ggsurvplots, c(list(x = psurvs, print = T), gridParams))
	dev.off()
	{% else %}
	for (ret in rets) {
		pvals   = rbind(pvals, ret$pvals)
		glbpv = rbind(glbpv, c(ret$var, ret$globalp))
		outpv   = file.path({{out.outdir | quote}}, paste0(ret$var, '.pval.txt'))
		outplot = file.path({{out.outdir | quote}}, paste0(ret$var, '.png'))
		write.table(ret$pvals, outpv, sep="\t", quote=F, col.names = T, row.names = F)
		do.call (png, c(outplot, devpars))
		print (ret$psurv)
		dev.off()
	}
	{% endif %}
	colnames(glbpv) = c('var', 'globalp')
	write.table(pvals, outpval, sep="\t", quote=F, col.names = T, row.names = F)
	write.table(glbpv[order(as.numeric(glbpv[,"globalp"])), , drop = F], glbpval, sep="\t", quote=F, col.names = T, row.names = F)

}
