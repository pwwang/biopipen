library(metap)
library(rlang)
library(dplyr)
library(data.table)
library(biopipen.utils)

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}
id_cols <- {{envs.id_cols | r}}
pval_col <- {{envs.pval_col | r}}
method <- {{envs.method | r}}
na <- {{envs.na | r}}
keep_single <- {{envs.keep_single | r}}
padj <- {{envs.padj | r}}

log <- get_logger()

if (method == "fisher") { method = "sumlog" }

# Check pval_cols
if (is.null(pval_col)) { stop("Must provide envs.pval_col") }

# Check id_cols
if (is.null(id_cols)) { stop("Must provide envs.id_cols") }
if (length(id_cols) == 1) {
    id_cols <- trimws(strsplit(id_cols, ",")[[1]])
}

log$info("Reading input and performing meta-analysis ...")
dt <- fread(infile, sep = "\t", header = TRUE, check.names = FALSE)
outdata <- dt[, list(N = .N, .pvals = list(get(pval_col))), by = id_cols]
outdata <- as.data.frame(outdata)

# Pre-allocate result vectors for maximum performance
n_rows <- nrow(outdata)
metaps <- numeric(n_rows)
ns <- integer(n_rows)

# Single-pass computation with pre-allocated vectors
for (i in seq_len(n_rows)) {
    ps <- outdata$.pvals[[i]]

    if (na == -1) {
        ps <- ps[!is.na(ps)]
    } else {
        ps[is.na(ps)] <- na
    }

    if (length(ps) == 0) {
        metaps[i] <- NA_real_
        ns[i] <- NA_integer_
    } else if (length(ps) == 1 && keep_single) {
        metaps[i] <- ps
        ns[i] <- 1L
    } else if (any(ps == 0)) {
        metaps[i] <- 0
        ns[i] <- length(ps)
    } else {
        metaps[i] <- do.call(method, list(ps))$p
        ns[i] <- length(ps)
    }
}

outdata$MetaPval <- metaps
outdata$N <- ns
outdata$.pvals <- NULL

# Sort by MetaPval
outdata <- setorder(as.data.table(outdata), MetaPval)
outdata <- as.data.frame(outdata)

if (padj != "none") {
    log$info("Calculating adjusted p-values ...")
    # Calculate adjusted p-values only on unique id_cols combinations
    dt <- as.data.table(outdata)
    pdata <- unique(dt, by = id_cols)[, c(id_cols, "MetaPval"), with = FALSE]
    pdata$MetaPadj <- p.adjust(pdata$MetaPval, method = padj)
    setkeyv(dt, id_cols)
    setkeyv(pdata, id_cols)
    outdata <- as.data.frame(dt[pdata[, c(id_cols, "MetaPadj"), with = FALSE]])
}

log$info("Writing output ...")
write.table(outdata, outfile, quote = FALSE, sep = "\t", row.names = FALSE)
