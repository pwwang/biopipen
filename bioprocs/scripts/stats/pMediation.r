library(methods)
library(mediation)
library(parallel)
{{rimport}}('__init__.r')
options(stringsAsFactors = FALSE)
set.seed(8525)

# inputs/outputs/arguments
infile   = {{i.infile | R}}
case0    = {{i.infile | fn2 | R}}
outfile  = {{o.outfile | R}}
outdir   = {{o.outdir | R}}
casefile = {{i.casefile | R}}
argscase = {{args.case | R}}
covfile  = {{args.cov | R}}
medopts  = {{args.medopts | R}}
inopts   = {{args.inopts | R}}
plotmed  = {{args.plot | R}}
pval     = {{args.pval | R}}
dofdr    = {{args.fdr | R}}
devpars  = {{args.devpars | R}}
nthread  = {{args.nthread | R}}

if (!is.list(plotmed) && plotmed == T) plotmed = list()
if (dofdr == T) dofdr = 'BH'

indata = read.table.inopts(infile, inopts)
if (casefile != "") {
	csdata = read.table(casefile, header = T, row.names = NULL, sep = "\t", check.names = F)
	colnames(csdata) = c('Case', 'Model', 'Fmula')
} else {
	if ('model' %in% names(argscase)) {
		argscase = list(argscase)
		names(argscase) = case0
	}
	csdata = as.data.frame(t(sapply(names(argscase), function(r) list(Case = r, Model = argscase[[r]]$model, Fmula = argscase[[r]]$fmula))))
}

covs = c()
if (covfile != '') {
	covdata = read.table(covfile, header = T, row.names = 1, sep = "\t", check.names = F)
	covs    = paste(bQuote(colnames(covdata)), collapse = '+')
	indata  = cbind(indata, covdata[rownames(indata)])
}

medanalysis = function(i) {
	tryCatch({
		# Case, Model, Fmula
		row = csdata[i,,drop = FALSE]
		logger('Dealing with case:', row$Case, '...')
		if (grepl(',', row$Model, fixed = TRUE)) {
			models = trimws(unlist(strsplit(row$Model, ',', fixed = T)))
			modelm = match.fun(models[2])
			modely = match.fun(models[1])
		} else {
			modelm = match.fun(row$Model)
			modely = match.fun(row$Model)
		}
		vars     = all.vars(as.formula(row$Fmula))
		outcome  = vars[1]
		treat    = vars[2]
		mediator = vars[3]
		fmulam   = as.formula(sprintf('%s ~ %s', bQuote(mediator), bQuote(treat)))
		fmulay   = as.formula(sprintf('%s ~ %s + %s', bQuote(outcome), bQuote(mediator), bQuote(treat)))
		if (length(covs) > 0) {
			fmulam = update.formula(fmulam, paste('. ~ +', covs))
			fmulay = update.formula(fmulay, paste('. ~ +', covs))
		}
		mm  = modelm(fmulam, data = indata)
		my  = modely(fmulay, data = indata)
		med = mediate(mm, my, treat = treat, mediator = mediator, outcome = outcome, boot = medopts$boot, sims = medopts$sims)
		if (is.na(med$d1.p) || is.na(med$n1)) {
			NULL
		} else {
			if (is.list(plotmed) && med$d1.p < pval && med$n1 > 0) {
				do.call(png, c(list(file.path(outdir, paste0(case, '.png'))), devpars))
				do.call(plot.mediate, c(list(med), plotmed))
				dev.off()
			}
			data.frame(
				Case      = row$Case,
				ACME      = med$d1,
				ACME95CI1 = med$d1.ci[1],
				ACME95CI2 = med$d1.ci[2],
				TotalE    = med$tau.coef,
				ADE       = med$z1,
				PropMed   = med$n1,
				Pval      = med$d1.p
			)
		}
	}, error = function(e) {
		logger(row$Case, ':', e, level = 'ERROR')
		NULL
	})
}

ret = do.call(rbind, mclapply(1:nrow(csdata), medanalysis, mc.cores = nthread))

if (is.null(ret)) {
	ret = data.frame(
		Case             = character(),
		ACME             = double(),
		ACME95CI1        = double(),
		ACME95CI2        = double(),
		TotalE           = double(),
		ADE              = double(),
		PropMed          = double(),
		Pval             = double()
	)
}

if (dofdr != F && nrow(ret) > 0) {
	ret$Qval = p.adjust(ret$Pval, method = dofdr)
}
if (nrow(ret)>0) {
	ret = ret[ret$Pval < pval & ret$PropMed > 0,,drop=F]
}

write.table(pretty.numbers(ret, list(ACME..ACME95CI1..ACME95CI2..PropMed..TotalE..ADE = '%.3f', Pval..Qval = '%.2E')), outfile, col.names = T, row.names = F, sep = "\t", quote = F)
