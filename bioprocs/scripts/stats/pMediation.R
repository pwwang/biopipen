{{"__init__.R" | rimport}}
library(methods)
library(mediation)
library(parallel)
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

if (!is.list(plotmed) && plotmed == TRUE) {
    plotmed = list()
}
if (dofdr == TRUE) {
    dofdr = 'BH'
}

indata = read.table.inopts(infile, inopts)
if (casefile != "") {
    csdata = read.table.inopts(casefile, list(rnames = TRUE, cnames = FALSE))
    if (ncol(csdata) == 1) {
        csdata$Model = "lm"
    } else {
        csdata = csdata[, 1:2, drop = FALSE]
    }
    colnames(csdata) = c('Fmula', 'Model')
} else {
    if ('model' %in% names(argscase)) {
        argscase = list(argscase)
        names(argscase) = case0
    }
    csdata = as.data.frame(t(sapply(
        names(argscase),
        function(r) list(Fmula = argscase[[r]]$fmula, Model = argscase[[r]]$model)
    )))
    rownames(csdata) = names(argscase)
}

covs = c()
if (covfile != '') {
    covdata = read.table.inopts(covfile, list(rnames = TRUE, cnames = TRUE))
    covs    = paste(bQuote(colnames(covdata)), collapse = '+')
}

medanalysis = function(case) {
    tryCatch({
        # Model, Fmula
        row = csdata[case,,drop = FALSE]
        logger('Dealing with case:', case, '...')
        modelm = match.fun(row$Model)
        modely = match.fun(row$Model)
        vars     = all.vars(as.formula(row$Fmula))
        outcome  = vars[1]
        treat    = vars[2]
        mediator = vars[3]
        fmulam   = as.formula(sprintf('%s ~ %s',
                                      bQuote(mediator),
                                      bQuote(treat)))
        fmulay   = as.formula(sprintf('%s ~ %s + %s',
                                      bQuote(outcome),
                                      bQuote(mediator),
                                      bQuote(treat)))
        mdata = indata
        if (length(covs) > 0) {
            fmulam = update.formula(fmulam, paste('. ~ +', covs))
            fmulay = update.formula(fmulay, paste('. ~ +', covs))
            mdata = cbind(mdata, covdata[rownames(mdata), ])
        }
        mm  = modelm(fmulam, data = mdata)
        my  = modely(fmulay, data = mdata)
        med = mediate(mm, my, treat = treat,
                      mediator = mediator,
                      outcome = outcome,
                      boot = medopts$boot,
                      sims = medopts$sims)
        if (is.na(med$d1.p) || is.na(med$n1)) {
            NULL
        } else {
            if (is.list(plotmed) && med$d1.p < pval && med$n1 > 0) {
                do.call(png, c(list(file.path(outdir, paste0(case, '.png'))), devpars))
                do.call(plot.mediate, c(list(med), plotmed))
                dev.off()
            }
            data.frame(
                Case      = case,
                Mediator  = mediator,
                Treat     = treat,
                Outcome   = outcome,
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
        logger(case, ':', e, level = 'ERROR')
        NULL
    })
}

ret = do.call(rbind, mclapply(rownames(csdata), medanalysis, mc.cores = nthread))

if (is.null(ret)) {
    ret = data.frame(
        Case      = character(),
        Mediator  = character(),
        Treat     = character(),
        Outcome   = character(),
        ACME      = double(),
        ACME95CI1 = double(),
        ACME95CI2 = double(),
        TotalE    = double(),
        ADE       = double(),
        PropMed   = double(),
        Pval      = double()
    )
}

if (dofdr != FALSE && nrow(ret) > 0) {
    ret$Qval = p.adjust(ret$Pval, method = dofdr)
}
if (nrow(ret)>0) {
    ret = ret[ret$Pval < pval & ret$PropMed > 0,,drop=F]
}

write.table(pretty.numbers(
    ret,
    list(ACME..ACME95CI1..ACME95CI2..PropMed..TotalE..ADE = '%.3f', Pval..Qval = '%.2E')
), outfile, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
