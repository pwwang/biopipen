library(rlang)
library(parallel)
library(mediation)
library(biopipen.utils)

infile <- {{in.infile | r}}
fmlfile <- {{in.fmlfile | r}}
outfile <- {{out.outfile | r}}

ncores <- {{envs.ncores | r}}
sims <- {{envs.sims | r}}
args <- {{envs.args | r}}
padj <- {{envs.padj | r}}
cases <- {{envs.cases | r}}
transpose_input <- {{envs.transpose_input | r}}

set.seed(123)
log <- get_logger()

log$info("Reading input file ...")
indata <- read.table(infile, header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE)
if (transpose_input) { indata <- t(indata) }

log$info("Reading formula file/cases ...")
if (!is.null(fmlfile)) {
    if (!is.null(cases) && length(cases) > 0) {
        log$warn("envs.cases ignored as in.fmlfile is provided")
    }
    fmldata <- read.table(fmlfile, header = TRUE, sep = "\t", row.names = NULL)
    # Case   M   Y   X   Cov     Model_M    Model_Y
    cases <- split(fmldata, fmldata$Case)
} else if (is.null(cases) || length(cases) == 0) {
    stop("Either envs.cases or in.fmlfile must be provided")
}

args <- args %||% list()

medanalysis <- function(i, total) {
    casename <- names(cases)[i]
    case <- cases[[casename]]
    if (total < 50) {
        log$info("- Case: ", casename)
    } else if (total < 500) {
        if (i %% 10 == 0) {
            log$info("- Processing case {i}/{total} ...")
        }
    } else {
        if (i %% 100 == 0) {
            log$info("- Processing case {i}/{total} ...")
        }
    }
    M <- case$M
    Y <- case$Y
    X <- case$X
    covs <- case$Cov
    modelm <- match.fun(case$Model_M)
    modely <- match.fun(case$Model_Y)
    fmlm <- as.formula(sprintf("%s ~ %s", bQuote(M), bQuote(X)))
    fmly <- as.formula(sprintf("%s ~ %s + %s", bQuote(Y), bQuote(M), bQuote(X)))
    if (!is.null(covs) && length(covs) == 1) {
        covs <- trimws(strsplit(covs, ",")[[1]])
    }
    if (!is.null(covs)) {
        cov_fml <- as.formula(sprintf("~ . + %s", paste(bQuote(covs), collapse = " + ")))
        fmlm <- update.formula(fmlm, cov_fml)
        fmly <- update.formula(fmly, cov_fml)
    }

    data <- indata[, c(M, X, Y, covs), drop = FALSE]
    data <- data[complete.cases(data), , drop = FALSE]
    margs <- args
    margs$sims <- sims
    margs$model.m <- modelm(fmlm, data = data)
    margs$model.y <- modely(fmly, data = data)
    margs$treat <- X
    margs$mediator <- M
    margs$outcome <- Y
    if (!is.null(covs)) {
        margs$covariates <- data[, covs, drop = FALSE]
    }
    med <- do_call(mediate, margs)
    if (is.na(med$d1.p) || is.na(med$n1)) {
        NULL
    } else {
        data.frame(
            Case         = casename,
            M            = M,
            X            = X,
            Y            = Y,
            ACME         = med$d1,
            ACME95CI1    = med$d1.ci[1],
            ACME95CI2    = med$d1.ci[2],
            TotalEffect  = med$tau.coef,
            ADE          = med$z1,
            PropMediated = med$n1,
            Pval         = med$d1.p
        )
    }
}

total <- length(cases)
out <- do_call(rbind, mclapply(1:total, medanalysis, total = total, mc.cores = ncores))

if (padj != "none") {
    out$Padj <- p.adjust(out$Pval, method = padj)
}

write.table(out, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
