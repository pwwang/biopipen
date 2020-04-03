{{"__init__.R" | rimport}}
library(parallel)
infile <- {{i.infile | quote}}
covfile <- {{i.covfile | quote}}
outfile <- {{o.outfile | quote}}
nthread <- {{args.nthread | R}}
tool <- {{args.tool | R}}
nk <- {{args.nk | R}}

#    g1 g2 g3
# s1 ...
# s2 ...
indata <- read.table.inopts(infile, list(rnames = TRUE, cnames = TRUE))
#    cov1 cov2
# s1 ...
# s2 ...
covdata <- read.table.inopts(covfile, list(rnames = TRUE, cnames = TRUE))
covdata <- covdata[colnames(indata), , drop = FALSE]
if (tool == 'peer') {
    library(peer)
    model <- PEER()
    PEER_setPhenoMean(model, as.matrix(indata))
    PEER_setNk(model, nk)
    PEER_setCovariates(model, as.matrix(covdata))
    PEER_setAdd_mean(model, TRUE)
    PEER_update(model)
    results <- PEER_getResiduals(model)
    colnames(results) <- colnames(indata)
    rownames(results) <- rownames(indata)
} else {
    indata <- t(indata)
    vars <- colnames(indata)
    covs <- colnames(covdata)

    regress_one <- function(var) {
        fmula <- as.formula(paste(bQuote(var), "~ ."))
        model <- lm(fmula, cbind(indata[, var, drop=FALSE], covdata))
        out <- t(as.data.frame(model$residuals))
        rownames(out) <- var
        out
    }

    results <- do.call(rbind, mcmapply(regress_one, vars, SIMPLIFY = FALSE, mc.cores = nthread))
}
write.table(results, outfile, sep = "\t", quote = FALSE)
