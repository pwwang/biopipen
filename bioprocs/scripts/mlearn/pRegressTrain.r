{{rimport}}('__init__.r')
options(stringsAsFactors = FALSE)

infile    = {{i.infile | quote}}
fmfile    = {{i.fmfile | quote}}
casefile  = {{i.casefile | quote}}
outfile   = {{o.outfile | quote}}
outdir    = {{o.outdir | quote}}
inopts    = {{args.inopts | R}}
args.plot = {{args.plot | R}}
covfile   = {{args.cov | R}}
devpars   = {{args.devpars | R}}
yval      = {{args.yval | R}}
glmfam    = {{args.glmfam}}
savemodel = {{args.save | R}}

regressone = function(idata, model, fmula, covs, case) {
	logger('Handing case:', case, '...')
	model = as.character(model)
	mdl   = eval(parse(text = model))
	ycol  = trimws(unlist(strsplit(fmula, "~", fixed = T))[1])
	ycol  = gsub('^`|`$', '', ycol)
	if (yval == 'categ') {
		ylabs   = levels(as.factor(idata[, ycol]))
		if (length(ylabs) < 2) 
			stop('Require more than 1 levels for categorical Y values in case ', case)
		ylevels = order(ylabs)
		ymin    = min(ylevels)
		ymax    = max(ylevels)
		idata[, ycol] = (match(idata[, ycol], ylabs) - ymin)/(ymax - ymin)
	} else {
		ymin = min(idata[, ycol], na.rm = T)
		ymax = max(idata[, ycol], na.rm = T)
		idata[, ycol] = (idata[, ycol] - ymin) / (ymax - ymin)
	}
	fml  = ifelse(
		is.null(covs), 
		as.character(fmula), 
		paste(c(as.character(fmula), sapply(covs, bQuote)), collapse = ' + ')
	)
	m = ifelse(
		model == 'glm', 
		mdl(as.formula(fml), data = idata, family = glmfam), 
		mdl(as.formula(fml), data = idata)
	)
	if (args.plot) {
		do.call(png, c(list(file.path(outdir, paste0(case, '.png'))), devpars))
		layout(matrix(1:4, ncol = 2))
		plot(m, pch = 20)
		dev.off()
	}
	m$ycol = ycol
	m$ymin = ymin
	m$ymax = ymax
	if (yval == 'categ') m$ylabs = ylabs
	if (savemodel)
		saveRDS(m, file.path(outdir, paste(case, model, 'rds', sep = '.')))
	ret = summary(m)$coefficients
	terms = rownames(ret)
	nterm = length(terms)
	cbind(Case = rep(case, nterm), Model = rep(model, nterm), Term = terms, ret)
}

indata = read.table.inopts(infile, inopts)
insts  = rownames(indata)
covs   = NULL
if (covfile!='') {
	covdata = read.table(covfile, header = T, row.names = 1, sep = "\t", check.names = F)
	covs    = colnames(covdata)
	indata  = cbind(indata, covdata[insts,,drop=F])
}

fmdata = read.table(fmfile, header = F, row.names = NULL, sep = "\t", check.names = F)
if (ncol(fmdata) < 3) 
	fmdata = cbind(Case = paste0('Case', 1:nrow(fmdata)), fmdata, stringsAsFactors = F)

casedata = NULL
if (casefile != '')
	casedata = read.table(casefile, header = F, row.names = NULL, sep = "\t", check.names = F)

cases = fmdata[,1]
out = data.frame(
	Case     = character(),
	Model    = character(),
	Term     = character(),
	Estimate = double(),
	StdError = double(),
	Stat     = double(),
	Pval     = double()
)
for (case in cases) {
	cdata = ifelse(
		is.null(casedata), 
		indata, 
		indata[casedata[which(casedata[,2]==case), 1],,drop=F])
	fdata = fmdata[which(fmdata[,1] == case),,drop=T]
	ret = regressone(cdata, fdata[[2]], fdata[[3]], covs, case)
	colnames(ret) = colnames(out)
	out = rbind(out, ret)
}
colnames(out) = c('Case', 'Model', 'Term', 'Estimate', 'StdError', 'Stat', 'Pval')

write.table(pretty.numbers(
	out, list(Estimate..StdError..Stat = '%.3f', Pval = '%.2E')
), outfile, row.names = F, col.names = T, sep = "\t", quote = F)
