library(methods)
{{rimport}}('__init__.r')
options(stringsAsFactors = F)

infile   = {{i.infile | R}}
gfile    = {{i.groupfile | R}}
cfile    = {{i.casefile | R}}
outfile  = {{o.outfile | R}}
outdir   = {{o.outdir | R}}
pcut     = {{args.pval | R}}
dofdr    = {{args.fdr | R}}
plotchow = {{args.plot | R}}
devpars  = {{args.devpars | R}}
ggs      = {{args.ggs | R}}
inopts   = {{args.inopts | R}}
covfile  = {{args.cov | R}}

if (plotchow) {
	{{rimport}}('plot.r')
}

if (dofdr == T) dofdr = 'BH'

chow.test = function(formula, group, data, covdata = NULL, ...) {
	fmvars  = all.vars(as.formula(formula))
	vars    = colnames(data)
	if (length(fmvars) == 2 && fmvars[2] == '.') {
		fmvars = c(fmvars[1], vars[vars!=fmvars[1] & vars!=group])
	}
	formula = sprintf(
		'%s ~ %s',
		bQuote(fmvars[1]),
		paste(sapply(fmvars[2:length(fmvars)], bQuote), collapse = '+')
	)
	covs = NULL
	if (is.null(covdata)) {
		pooledfm = as.formula(formula)
	} else {
		covdata = covdata[rownames(data),,drop = FALSE]
		covs    = colnames(covdata)
		data    = cbind(data, covdata)
		rm(covdata)
		pooledfm = as.formula(paste(formula, '+', paste(sapply(covs, bQuote), collapse = '+')))
		fmvars   = c(fmvars, covs)
	}
	
	if (sum(complete.cases(data[,fmvars])) < 2) {
		pooled_lm = NULL
	} else {
		pooled_lm = lm(pooledfm, data = data, ...)
	}
	#coeff     = as.list(pooled_lm$coefficients)
	groups    = levels(as.factor(data[, group]))
	group_lms = sapply(groups, function(g) {
		if (is.null(covdata)) {
			subfm  = as.formula(formula)
		} else {
			subfm = as.formula(paste(
				formula, '+', 
				#paste(sapply(covs, function(x) paste0('offset(', coeff[[x]], '*', bQuote(x), ')')), collapse = '+')
				paste(sapply(covs, bQuote), collapse = '+')
			))
		}
		sublmdata = data[data[,group] == g, , drop = FALSE]
		if (sum(complete.cases(sublmdata[,fmvars])) < 2) {
			NULL
		} else {
			list(lm(subfm, data = sublmdata, ...))
		}
	})

	pooled.ssr = ifelse(is.null(pooled_lm), NA, sum(pooled_lm$residuals ^ 2))
	subssr     = ifelse(is.false(group_lms, 'any'), NA, sum(sapply(group_lms, function(x) sum(x$residuals ^ 2))))
	ngroups    = length(groups)
	K          = length(fmvars) + length(covs)
	J          = (ngroups - 1) * K
	DF         = nrow(data) - ngroups * K
	FS         = (pooled.ssr - subssr) * DF / subssr / J
	list(
		pooled.lm  = pooled_lm,
		group.lms  = group_lms,
		Fstat      = FS,
		group      = group,
		pooled.ssr = pooled.ssr,
		group.ssr  = subssr,
		Pval       = pf(FS, J, DF, lower.tail = FALSE)
	)
}

plot.chow = function(chow, plotfile, ggs, devpars) {
	cols     = all.vars(chow$pooled.lm$terms)[1:2]
	plotdata = do.call(rbind, lapply(names(chow$group.lms), function(m) data.frame(chow$group.lms[[m]]$model[, cols, drop = FALSE], group = m)))
	colnames(plotdata)[3] = chow$group
	if (!is.null(ggs$scale_color_discrete)) {
		ggs$scale_color_discrete$name = ifelse(
			is.function(ggs$scale_color_discrete$name), 
			ggs$scale_color_discrete$name(chow$group),
			chow$group
		)
		ggs$scale_color_discrete$labels = sapply(names(chow$group.lms), function(x) {
			coeff = as.list(chow$group.lms[[x]]$coefficients)
			bquote(.(x): beta[.(cols[2])]==.(round(coeff[[cols[2]]], 3)) ~ "," ~ epsilon == .(round(coeff[['(Intercept)']], 3)))
		})
	}

	if (!is.null(ggs$scale_shape_discrete)) {
		ggs$scale_shape_discrete$name = ifelse(
			is.function(ggs$scale_shape_discrete$name), 
			ggs$scale_shape_discrete$name(chow$group),
			chow$group
		)
		ggs$scale_shape_discrete$labels = sapply(names(chow$group.lms), function(x) {
			coeff = as.list(chow$group.lms[[x]]$coefficients)
			bquote(.(x): beta[.(cols[2])]==.(round(coeff[[cols[2]]], 3)) ~ "," ~ epsilon == .(round(coeff[['(Intercept)']], 3)))
		})
	}

	if (is.null(ggs$scale_color_discrete) && is.null(ggs$scale_shape_discrete)) {
		ggs$scale_color_discrete = list(
			name = chow$group,
			labels = sapply(names(chow$group.lms), function(x) {
				coeff = as.list(chow$group.lms[[x]]$coefficients)
				bquote(.(x): beta[.(cols[2])]==.(round(coeff[[cols[2]]], 3)) ~ "," ~ epsilon == .(round(coeff[['(Intercept)']], 3)))
			})
		)
	}

	if (is.null(ggs$scale_color_discrete)) ggs$scale_color_discrete = ggs$scale_shape_discrete
	if (is.null(ggs$scale_shape_discrete)) ggs$scale_shape_discrete = ggs$scale_color_discrete
	
	plot.points(
		plotdata, 
		plotfile, 
		x = 2, y = 1, 
		params = list(aes_string(color = chow$group, shape = chow$group)), 
		ggs = c(ggs, list(
			geom_smooth = list(aes_string(color = chow$group), method = "lm", se = FALSE)
		))
	)
}

formatlm = function(m) {
	if (class(m) == 'lm') {
		coeff = as.list(m$coefficients)
		vars = all.vars(m$terms)
		terms = unlist(sapply(c(vars[2:length(vars)], '(Intercept)', 'N'), function(x) {
			ce = list.get(coeff, x, list.get(coeff, bQuote(x)))
			if (x == 'N') {
				paste0('N=', nrow(m$model))
			} else if (is.null(ce)) {
				NULL
			} else {
				l = ifelse(x == '(Intercept)', '_', x)
				paste0(l, '=', round(ce, 3))
			}
		}))
		paste(terms[!is.null(terms)], collapse = ', ')
	} else {
		paste(sapply(names(m), function(x) {
			paste0(x, ': ', formatlm(m[[x]]))
		}), collapse = ' // ')
	}
}

results = data.frame(
	Case   = character(),
	Pooled = character(),
	Groups = character(),
	SSR    = double(),
	SumSSR = double(),
	Fstat  = double(),
	Pval   = double()
)

indata = read.table.inopts(infile, inopts, try = TRUE)
if (is.null(indata)) {
	write.table(results, outfile, col.names = T, row.names = F, sep = "\t", quote = F)
	quit(save = "no")
}

#     X1  X2  X3  X4 ... Y
# G1  1   2   1   4  ... 9
# G2  2   3   1   1  ... 3
# ... ...
# Gm  3   9   1   7  ... 8
#K      = ncol(indata)
covdata = NULL
covs    = NULL
if (covfile != "") {
	covdata = read.table(covfile, header = T, row.names = 1, check.names = F)
	#indata  = cbind(covdata[rownames(indata),,drop = F], indata)
	covs = colnames(covs)
}
gdata  = read.table.inopts(gfile, list(cnames = TRUE, rnames = TRUE))
# 	Case1	Case2
# G1	Group1	Group1
# G2	Group1	NA
# G3	Group2	Group1
# ... ...
# Gm	Group2	Group2
cases  = colnames(gdata)
fmulas = data.frame(x = rep(paste(bQuote(colnames(indata)[ncol(indata)]), '~ .'), length(cases)))
rownames(fmulas) = cases
if (!is.null(cfile) && cfile != "") {
	fmulas = read.table.inopts(cfile, list(cnames = FALSE, rnames = TRUE))
	cases  = rownames(fmulas)
}

for (case in cases) {
	logger('Handling case: ', case, '...')
	fmula  = fmulas[case,,drop = TRUE]
	groups = gdata[!is.na(gdata[,case]),case,drop = FALSE]
	data   = cbind(indata[rownames(groups),, drop = FALSE], group = groups)
	colnames(data)[ncol(data)] = case
	ct = chow.test(fmula, case, data, covdata = covdata)
	if (dofdr == FALSE && (is.na(ct$Pval) || ct$Pval >= pcut)) {
		next
	}
	results = rbind(results, list(
		Case   = case,
		Pooled = formatlm(ct$pooled.lm),
		Groups = formatlm(ct$group.lms),
		SSR    = ct$pooled.ssr,
		SumSSR = ct$group.ssr,
		Fstat  = ct$Fstat,
		Pval   = ct$Pval
	))
	# doplot
	if (plotchow && ct$Pval < pcut) {
		plot.chow(ct, file.path(outdir, paste0(case, '.png')), ggs, devpars)
	}
}

if (dofdr != F) {
	results = cbind(results, Qval = p.adjust(results$Pval, method = dofdr))
} 
write.table(pretty.numbers(results, list(
	SSR..SumSSR..Fstat = '%.3f',
	Pval..Qval = '%.3E'
)), outfile, col.names = T, row.names = F, sep = "\t", quote = F)

