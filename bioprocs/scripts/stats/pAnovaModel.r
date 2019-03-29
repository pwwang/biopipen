{{rimport}}('__init__.r')

infile   = {{i.infile | quote}}
casefile = {{i.casefile | quote}}
outfile  = {{o.outfile | quote}}
model    = {{args.model | quote}}
fmula    = {{args.fmula | R}}
inopts   = {{args.inopts | R}}
cov      = {{args.cov | quote}}

indata  = read.table.inopts(infile, inopts)
covs    = NULL
if (cov != '') {
	covdata = read.table.inopts(cov, list(cnames = TRUE, rnames = TRUE))
	covs    = colnames(covdata)
	indata  = cbind(indata, covdata[rownames(indata),,drop = FALSE])
} else if (is.null(fmula)) {
	stop('Either `i.casefile` or `args.fmula` should be specified.')
}

add.covs = function(fmula) {
	if (is.null(covs) || is.null(fmula)) {
		return(fmula)
	} else {
		asfmula = as.formula(fmula)
		allvars = all.vars(asfmula)
		if (length(allvars) == 2 && allvars[2] == '.') {
			return (fmula)
		} else {
			return (paste(fmula, paste(sapply(covs, bQuote), collapse = '+'), sep = '+'))
		}
	}
}

do.case = function(case, m, fmula1, fmula2 = NULL) {
	fmula1 = as.formula(add.covs(fmula1))
	fmula2 = ifelse(is.null(fmula2), NULL, as.formula(add.covs(fmula2)))
	model1 = do.call(m, list(fmula1, data = indata))
	if (!is.null(fmula2)) {
		model2 = do.call(m, list(fmula2, data = indata))
		ret = anova(model2, mode1)
		colnames(ret) = c('Res_DF', 'RSS', 'DF', 'SumSq', 'Fstat', 'Pval')
		ret = cbind(list(Case = case), ret[2,,drop = FALSE])
	} else {
		ret = anova(model1)
		colnames(ret) = c('DF', 'SumSq', 'MeanSq', 'Fstat', 'Pval')
		ret = ret[complete.cases(ret),,drop = FALSE]
		ret = cbind(list(Case = case), list(Var = rownames(ret)), ret)
		rownames(ret) = NULL
	}
	ret
}

if (casefile == '') {
	casedata = data.frame(Case = 'Case1', Formula = fmula[1], Formula2 = fmula[2], Model = model)
} else {
	casedata = read.table.inopts(casefile, list(rnames = FALSE, cnames = TRUE))
	cdncol = ncol(casedata)
	if (cdncol < 2) {
		stop('Case file requires >=2 columns (Case and Formula)')
	} else if (cdncol == 2) {
		colnames(casedata) = c('Case', 'Formula')
	} else if (cdncol == 3) {
		colnames(casedata) = c('Case', 'Formula', 'Formula2')
	} else {
		colnames(casedata) = c('Case', 'Formula', 'Formula2', 'Model')
	}
	if (!'Formula2' %in% colnames(casedata)) {
		casedata$Formula2 = NA
	}
	if (!'Model' %in% colnames(casedata)) {
		casedata$Model = 'lm'
	}
}

r = do.call(rbind, apply(casedata, 1, function(row) {
	do.case(row[1], row[4], row[2], ifelse(is.na(row[3]), NULL, row[3]))
}))

write.table(r, outfile, quote = FALSE, sep = "\t")

