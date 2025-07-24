library(dcanr)
library(scuttle)
library(doRNG)
library(doParallel)
library(snpStats)
library(rlang)
library(dplyr)
library(biopipen.utils)

infile <- {{in.infile | r}}
groupfile <- {{in.groupfile | r}}
outfile <- {{out.outfile | r}}
method <- {{envs.method | r}}
beta <- {{envs.beta | r}}
padj <- {{envs.padj | r}}
perm_batch <- {{envs.perm_batch | r}}
seed <- {{envs.seed | r}}
ncores <- {{envs.ncores | r}}
transpose_input <- {{envs.transpose_input | r}}
transpose_group <- {{envs.transpose_group | r}}

log <- get_logger()

log$info("Setting seed and parallel backend ...")
set.seed(seed)
registerDoParallel(cores = ncores)
registerDoRNG(seed)

log$info("Reading input files ...")
indata <- read.table(infile, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
if (transpose_input) {
    indata <- t(indata)
}
gdata <- read.table(groupfile, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
if (transpose_group) {
    gdata <- t(gdata)
}
ngroups <- ncol(gdata)

sign2 <- function(x) sign(x) * x^2
mat2vec <- dcanr:::mat2vec

diffcoex_score <- function(group) {

    gvals <- unique(gdata[, group, drop = TRUE])
    if (length(gvals) < 2) {
        log$debug("  Less than 2 groups in the input. Skipping ...")
        return(NULL)
    }
    rs <- lapply(gvals, function(gval) {
        samples <- rownames(gdata[gdata[[group]] == gval, , drop = FALSE])
        expr <- indata[samples, , drop = FALSE]
        if (length(samples) < 3) {
            log$debug("  Less than 3 samples in one of the groups. Skipping ...")
            return(NULL)
        }
        cor.pairs(as.matrix(expr), cor.method = method)
    })
    rs[sapply(rs, is.null)] <- NULL
    if (length(rs) < 2) {
        log$debug("  Less than 2 groups with at least 3 samples. Skipping ...")
        return(NULL)
    }
    N <- length(rs)
    C0 <- lapply(rs, sign2)
    C0 <- Reduce(`+`, C0) / N
    D <- lapply(rs, function(r) abs(sign2(r) - C0))
    D <- Reduce(`+`, D) / 2 / (N - 1)
    D <- sqrt(D)
    D <- D^beta
    T_ovlap <- D %*% D + ncol(D) * D  #calc topological ovlap

    mins = matrix(rep(rowSums(D), ncol(D)), nrow = ncol(D))
    mins = pmin(mins, matrix(rep(colSums(D), each = ncol(D)), nrow = ncol(D)))
    T_ovlap = 1 - (T_ovlap/(mins + 1 - D))

    diag(T_ovlap) = 1

    #add run parameters as attributes
    attributes(T_ovlap) = c(
        attributes(T_ovlap),
        'method' = method,
        'beta' = beta,
        'call' = match.call()
    )

    return(1 - T_ovlap)
}


perm_test <- function(dcscores, group, B = perm_batch) {
    obs = mat2vec(dcscores)

    #package requirements
    pckgs = c('dcanr')

    #perform permutation
    pvals = foreach(
        b = seq_len(B),
        .combine = function(...) {mapply(sum, ...)},
        .multicombine = TRUE,
        .inorder = FALSE,
        .packages = pckgs
    ) %dorng% {
        #shuffle condition and recalculate scores
        env = new.env()
        assign('group', group, envir = env)
        permsc = eval(attr(dcscores, 'call'), envir = env)
        permsc = mat2vec(permsc)

        #count elements greater than obs
        permsc = abs(permsc)
        permsc = permsc[!(is.na(permsc) || is.infinite(permsc))]
        permcounts = vapply(abs(obs), function(x) sum(permsc > x), 0)
        return(c(permcounts, length(permsc)))
    }

    #p-values
    N <- pvals[length(pvals)]
    pvals <- pvals[-(length(pvals))] / N
    # attributes(pvals) = attributes(obs)
    # pvals = dcanr:::vec2mat(pvals)
    # attr(pvals, 'dc.test') = 'permutation'
    # return(pvals)
    # Format into Group,Feature1,Feature2,Pval
    feature_pairs <- as.data.frame(t(combn(attr(obs, 'feature.names'), 2)))
    colnames(feature_pairs) <- c('Feature1', 'Feature2')
    feature_pairs$Group <- group
    feature_pairs$Pval <- pvals
    feature_pairs[, c('Group', 'Feature1', 'Feature2', 'Pval'), drop = FALSE]
}

do_one_group <- function(i) {
    group <- colnames(gdata)[i]
    log$info("- Processing group {i}/{ngroups}: {group} ...")
    log$info("  Calculating differential co-expression scores ...")
    dcscores <- diffcoex_score(group)

    if (!is.null(dcscores)) {
        log$info("  Calculating p-values ...")
        perm_test(dcscores, group)
    }
}

trios <- do_call(rbind, lapply(seq_len(ngroups), do_one_group))
if (padj != "none") {
    log$info("Correcting p-values ...")
    trios$Padj <- p.adjust(trios$Pval, method = padj)
}

log$info("Writing output ...")
write.table(trios, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
