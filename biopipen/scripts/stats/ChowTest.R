source("{{biopipen_dir}}/utils/misc.R")

library(rlang)
library(dplyr)

infile <- {{in.infile | r}}
groupfile <- {{in.groupfile | r}}
fmlfile <- {{in.fmlfile | r}}
outfile <- {{out.outfile | r}}
padj <- {{envs.padj | r}}
transpose_input <- {{envs.transpose_input | r}}
transpose_group <- {{envs.transpose_group | r}}

log_info("Reading input files ...")
indata <- read.table(infile, header = TRUE, sep = "\t", row.names = 1)
if (transpose_input) {
	indata <- t(indata)
}
groupdata <- read.table(groupfile, header = TRUE, sep = "\t", row.names = 1)
if (transpose_group) {
	groupdata <- t(groupdata)
}
fmldata <- read.table(fmlfile, header = TRUE, sep = "\t", row.names = NULL)
colnames(fmldata)[1:2] <- c("Group", "Formula")

chow.test <- function(fml, grouping) {
    formula <- as.formula(fml)
    pooled_lm <- tryCatch(lm(formula, data = indata), error = function(e) NULL)
    if (is.null(pooled_lm)) {
        return(list(
            pooled.lm  = NA,
            group.lms  = NULL,
            Fstat      = NA,
            group      = grouping,
            pooled.ssr = NA,
            group.ssr  = NA,
            Pval       = NA
        ))
    }

    splitdata <- split(indata, groupdata[rownames(indata), grouping])
    group_lms <- lapply(names(splitdata), function(g) {
        tryCatch(lm(formula, data = splitdata[[g]]), error = function(e) NULL)
    })
	names(group_lms) <- names(splitdata)

    fmvars     <- all.vars(formula)
	pooled.ssr <- sum(pooled_lm$residuals ^ 2)
	subssr     <- ifelse(any(is.null(group_lms)), NA, sum(sapply(group_lms, function(x) sum(x$residuals ^ 2))))
	ngroups    <- length(splitdata)
	K          <- ifelse(fmvars[2] == ".", ncol(indata), length(fmvars))
	J          <- (ngroups - 1) * K
	DF         <- nrow(indata) - ngroups * K
	FS         <- (pooled.ssr - subssr) * DF / subssr / J
	list(
		pooled.lm  = pooled_lm,
		group.lms  = group_lms,
		Fstat      = FS,
		group      = grouping,
		pooled.ssr = pooled.ssr,
		group.ssr  = subssr,
		Pval       = pf(FS, J, DF, lower.tail = FALSE)
	)
}

formatlm <- function(m) {
	if (class(m) == 'lm') {
		coeff <- as.list(m$coefficients)
		vars <- all.vars(m$terms)
		terms <- unlist(sapply(na.omit(c(vars[2:length(vars)], '(Intercept)', 'N')), function(x) {
			ce <- coeff[[x]] %||% coeff[[bQuote(x)]]
			if (x == 'N') {
				paste0('N=', nrow(m$model))
			} else if (is.null(ce)) {
				NULL
			} else {
				l <- ifelse(x == '(Intercept)', '_', x)
				paste0(l, '=', round(ce, 3))
			}
		}))
		paste(terms[!is.null(terms)], collapse = ', ')
	} else {
		paste(sapply(names(m), function(x) {
			paste0(x, ': ', formatlm(m[[x]]))
		}), collapse = ' // ')
	}
}

log_info("Running Chow tests ...")
ncases <- nrow(fmldata)
results <- do_call(rbind, lapply(
    seq_len(ncases),
    function(i) {
		fmlrow <- fmldata[i, , drop=TRUE]
        if (i %% 100 == 0) {
            log_info("- {i} / {ncases} ...")
        }
        log_debug("  Running Chow test for formula: {fmlrow$Formula} (grouping = {fmlrow$Group})")

        res <- chow.test(fmlrow$Formula, fmlrow$Group)
		fmlrow$Pooled <- formatlm(res$pooled.lm)
		fmlrow$Groups <- formatlm(res$group.lms)
		fmlrow$SSR <- res$group.ssr
		fmlrow$SumSSR <- res$pooled.ssr
		fmlrow$Fstat <- res$Fstat
		fmlrow$Pval <- res$Pval
		fmlrow
    }
)) %>% as.data.frame()

if (padj != "none") {
    log_info("Adjusting p-values ...")
    results$Padj <- p.adjust(results$Pval, method = padj)
}

log_info("Writing output ...")
# unimplemented type 'list' in 'EncodeElement'
results <- apply(results, 2, as.character)
write.table(results, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
