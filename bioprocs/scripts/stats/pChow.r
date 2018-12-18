library(methods)
{{rimport}}('plot.r', '__init__.r')
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
if (dofdr == T) dofdr = 'BH'

regress = function(regdata, name, fmula = NULL) {
	Y = bQuote(colnames(regdata)[ncol(regdata)])
	fmula = ifelse(is.null(fmula), paste(Y, "~ ."), fmula)
	m = lm(as.formula(fmula), data = regdata)
	list(model = m, ssr = sum(m$residuals ^ 2), n = nrow(regdata), name = name)
}

formlm = function(model, k, withname = T) {
	lencoef = length(model$model$coefficients)
	coefns  = c(names(model$model$coefficients)[(lencoef-k+2):lencoef], '_')
	coeffs  = as.numeric(c(model$model$coefficients[(lencoef-k+2):lencoef], model$model$coefficients[1]))
	if (withname) {
		paste0(model$name, ': ', paste(coefns, round(coeffs, 3), sep = '=', collapse = ' '))
	} else {
		paste(coefns, round(coeffs, 3), sep = '=', collapse = ' ')
	}
}

chow = function(pooled, subregs, k, case) {
	subssr  = sum(sapply(subregs, function(x) x$ssr))
	ngroups = length(subregs)
	J  = (ngroups - 1) * k
	DF = pooled$n - ngroups * k
	FS = (pooled$ssr - subssr) * DF / subssr / J
	groups = lapply(subregs, function(m) formlm(m, k))
	pooledm = formlm(pooled, k, FALSE)
	list(Case = case, Pooled = pooledm, Groups = paste(groups, collapse = '; '), Fstat = FS, Pval = pf(FS, J, DF, lower.tail = FALSE))
}

model2eq = function(model) {
	vars = colnames(model$model)
	cf   = sapply(model$coefficients, function(f) {
		if (is.na(f)) return("NA")
		f = round(f, 2)
		if (f >= 0) return(paste0("+", f))
		return(paste0("-", -f))
	})
	paste0(vars[1], ' = ', cf[1], paste(
		sapply(
			2:length(cf),
			function(x) paste0(cf[x], '*', vars[x])
		), collapse = ''
	))
}

indata = read.table.inopts(infile, inopts)
#     X1  X2  X3  X4 ... Y
# G1  1   2   1   4  ... 9
# G2  2   3   1   1  ... 3
# ... ...
# Gm  3   9   1   7  ... 8
#K = ncol(indata)
if (covfile != "") {
	covdata = read.table(covfile, header = T, row.names = 1, check.names = F)
	indata  = cbind(covdata[rownames(indata),,drop = F], indata)
}
gdata  = read.table(gfile, header = F, row.names = NULL, check.names = F)
# G1	Group1	Case1
# G2	Group1	Case1
# ... ...
# Gs	Group2	Case1
# Gt	Group2	Case1
# Gt	Group1	Case2
# ... ...
# Gu	Group1	Case2
# ... ...
# Gz	Group2	Case2
if (ncol(gdata) == 2) { # no case
	cases = 'Case1'
}
fmulas = NULL
if (cfile != "") {
	fmulas = read.table(cfile, header = F, row.names = 1, sep = "\t", check.names = F)
	cases  = rownames(fmulas)
}

results = data.frame(
	Case   = character(),
	Pooled = character(),
	Groups = character(),
	Fstat  = double(),
	Pval   = double())
for (case in cases) {
	caserows  = if(ncol(gdata) > 2) gdata[which(gdata$Case == case),,drop = F] else gdata
	fmula     = fmulas[case,]
	K         = ifelse(is.null(fmula), ncol(indata), length(unlist(strsplit(fmula, "\\s*(\\~|\\+)\\s*"))))
	pooled_lm = regress(indata[caserows[,1],,drop = F], name = 'Pooled', fmula = fmulas[case,])
	subgroups = levels(factor(caserows[,2]))
	subgrp_lm = lapply(subgroups, function(sgroup) {
		sdata = indata[caserows[which(caserows[,2] == sgroup),1],,drop = F]
		if (nrow(sdata) < 3) NULL
		else regress(sdata, name = sgroup, fmula = fmulas[case,])
	})
	subgrp_lm[sapply(subgrp_lm, is.null)] <- NULL
	# no subgroups
	if (length(subgrp_lm) < 2) next
	ret = chow(pooled_lm, subgrp_lm, k = K, case = case)
	if (is.na(ret$Pval) || ret$Pval >= pcut) next
	results = rbind(results, ret)

	# doplot
	if (plotchow) {
		incol = ncol(indata)
		plotdata = cbind(indata[caserows[,1],,drop = F], Group = caserows[,2])
		rcase = make.names(case)
		colnames(plotdata)[incol + 1] = rcase

		labels = sapply(c(subgrp_lm, list(pooled_lm)), function(m) paste0(m$name, ': ', model2eq(m$model)))
		ggs1 = c(ggs, list(
			geom_smooth = list(
				aes(color = 'Pooled'),
				method   = 'lm',
				se       = F,
				linetype = "twodash"
			),
			geom_smooth = list(
				aes_string(color = rcase),
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
				values = c(scales::hue_pal()(length(subgroups)), '#555555'), 
				name   = "",
				limit  = c(subgroups, 'Pooled'),
				labels = labels
			)
		))
		plot.scatter(
			plotdata, 
			file.path(outdir, paste0(case, '.png')), 
			x      = incol - 1,
			y      = incol,
			ggs    = ggs1,
			params = list(aes_string(shape = rcase, color = rcase))
		)
	}
}

if (dofdr != F) {
	results = cbind(results, Qval = p.adjust(results$Pval, method = dofdr))
} 
write.table(pretty.numbers(results, list(
	Fstat      = '%.3f',
	Pval..Qval = '%.3E'
)), outfile, col.names = T, row.names = F, sep = "\t", quote = F)

