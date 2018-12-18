library(methods)
library(mediation)
{{rimport}}('__init__.r')

set.seed(8525)

infile   = {{i.infile | R}}
outfile  = {{o.outfile | R}}
outdir   = {{o.outdir | R}}
medfile  = {{i.medfile | R}}
casefile = {{i.casefile | R}}
covfile  = {{args.cov | R}}
medopts  = {{args.medopts | R}}
inopts   = {{args.inopts | R}}
plotmed  = {{args.plot | R}}
pval     = {{args.pval | R}}
dofdr    = {{args.fdr | R}}
devpars  = {{args.devpars | R}}

if (!is.list(plotmed) && plotmed == T) plotmed = list()
if (dofdr == T) dofdr = 'BH'

getfmula = function(Y, feats) {
	sprintf("%s ~ %s", bQuote(Y), paste(sapply(feats, bQuote), collapse = "+"))
}

medanalysis = function(idata, mediator, treat, outcome, modelm, modely, fmulam, fmulay) {
	covs   = setdiff(colnames(idata), c(mediator, treat, outcome))
	modelm = eval(parse(text = modelm))
	modely = eval(parse(text = modely))
	fmulam = paste(c(fmulam, covs), collapse = ' + ')
	fmulay = paste(c(fmulay, covs), collapse = ' + ')
	mm     = modelm(as.formula(fmulam), data = idata)
	my     = modely(as.formula(fmulay), data = idata)
	mediate(mm, my, treat = treat, mediator = mediator, outcome = outcome, boot = medopts$boot, sims = medopts$sims)
}

indata = read.table.inopts(infile, inopts)
# Case	Mediator	Treat	Outcome	MedelM	ModelY	FmulaM	FmulaY
medata = read.table(medfile, header = T, row.names = NULL, sep = "\t", check.names = F, stringsAsFactors = F)

medcols = colnames(medata)
insts   = rownames(indata)
#if ('Case' %in% medcols && casefile == '') {
#	stop('No casefile provided with "Case" specified in medfile.')
#}

if (covfile != '') {
	covdata = read.table(covfile, header = T, row.names = 1, sep = "\t", check.names = F)
	indata  = cbind(indata, covdata[rownames(indata)])
}

specases = NULL
if (casefile != '') {
	specases = read.table(casefile, header = F, row.names = NULL, check.names = F, sep = "\t")
}

rownames(medata) = paste0('Case', 1:nrow(medata))
if ('Case' %in% medcols) {
	rownames(medata) = medata$Case
	medata = medata[, -which('Case' %in% medcols), drop = F]
}

ret = data.frame(
	Case             = character(),
	ACME             = double(),
	ACME95CI1        = double(),
	ACME95CI2        = double(),
	TotalE           = double(),
	ADE              = double(),
	PropMed          = double(),
	Pval             = double(),
	stringsAsFactors = F)
for (case in rownames(medata)) {
	rows    = if(is.null(specases)) insts else specases[which(specases[,2] == case), 1]
	medinfo = medata[case,,drop = F]
	med     = medanalysis(indata[rows,, drop = F], medinfo$Mediator, medinfo$Treat, medinfo$Outcome, medinfo$ModelM, medinfo$ModelY, medinfo$FmulaM, medinfo$FmulaY)
	ret = rbind(ret, list(
		Case      = case,
		ACME      = med$d1,
		ACME95CI1 = med$d1.ci[1],
		ACME95CI2 = med$d1.ci[2],
		TotalE    = med$tau.coef,
		ADE       = med$z1,
		PropMed   = med$n1,
		Pval      = med$d1.p
	), stringsAsFactors = F)
	if (is.na(med$d1.p) || med$d1.p >= pval || is.na(med$n1) || med$n1 <= 0) next
	if (is.list(plotmed)) {
		do.call(png, c(list(file.path(outdir, paste0(case, '.png'))), devpars))
		do.call(plot.mediate, c(list(med), plotmed))
		dev.off()
	}
}
if (dofdr != F && nrow(ret) > 0) {
	ret$Qval = p.adjust(ret$Pval, method = dofdr)
}

write.table(pretty.numbers(ret, list(ACME..ACME95CI1..ACME95CI2..PropMed..TotalE..ADE = '%.3f', Pval..Qval = '%.2E')), outfile, col.names = T, row.names = F, sep = "\t", quote = F)
