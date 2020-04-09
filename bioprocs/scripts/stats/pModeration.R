{{"__init__.R" | rimport}}
library(methods)
library(parallel)
library(rockchalk)
options(stringsAsFactors = FALSE)
set.seed(8525)

# inputs/outputs/arguments
infile   = {{i.infile | R}}
case0    = {{i.infile | fn2 | R}}
outfile  = {{o.outfile | R}}
outdir   = {{o.outdir | R}}
casefile = {{i.casefile | R}}
argscase = {{args.case | R}}
inopts   = {{args.inopts | R}}
plotmod  = {{args.plot | R}}
pval     = {{args.pval | R}}
dofdr    = {{args.fdr | R}}
devpars  = {{args.devpars | R}}
nthread  = {{args.nthread | R}}

if (!is.list(plotmod) && plotmod == TRUE) {
    plotmod = list()
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

modanalysis = function(case) {
    tryCatch({
        # Model, Fmula
        row = csdata[case,,drop = FALSE]
        logger('Dealing with case:', case, '...')
        model = match.fun(row$Model)
        vars = all.vars(as.formula(row$Fmula))
        outcome = vars[1]
        treat = vars[2]
        moderator = vars[3]
        fmula = as.formula(sprintf('%s ~ %s*%s',
                                   bQuote(outcome),
                                   bQuote(moderator),
                                   bQuote(treat)))
        mm = model(fmula, data = indata)
        mod = as.data.frame(summary(mm)$coefficients[4,,drop=FALSE])
        colnames(mod) = c("Coeff", "StdErr", "Tval", "Pval")

        if (is.na(mod$Pval)) {
            NULL
        } else {
            if (is.list(plotmod) && mod$Pval < pval) {
                do.call(png, c(list(
                    file.path(outdir,
                              paste0(case, '-', moderator, '-', treat, '-', outcome, '.png'))
                ), devpars))
                plotSlopes(mm, plotx = treat, modx = moderator, xlab = treat, ylab = outcome)
                dev.off()
            }
            data.frame(
                Case      = case,
                Moderator = moderator,
                Treat     = treat,
                Outcome   = outcome,
                Coeff     = mod$Coeff,
                StdErr    = mod$StdErr,
                Tval      = mod$Tval,
                Pval      = mod$Pval
            )
        }
    }, error = function(e) {
        logger(case, ':', e, level = 'ERROR')
        NULL
    })
}

ret = do.call(rbind, mclapply(rownames(csdata), modanalysis, mc.cores = nthread))

if (is.null(ret)) {
    ret = data.frame(
        Case      = character(),
        Moderator = character(),
        Treat     = character(),
        Outcome   = character(),
        Coeff     = double(),
        StdErr    = double(),
        Tval      = double(),
        Pval      = double()
    )
}

if (dofdr != FALSE && nrow(ret) > 0) {
    ret$Qval = p.adjust(ret$Pval, method = dofdr)
}
if (nrow(ret)>0) {
    ret = ret[ret$Pval < pval,,drop=F]
}

write.table(pretty.numbers(
    ret,
    list(Coeff..StdErr..Tval = '%.3f', Pval = '%.2E')
), outfile, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
