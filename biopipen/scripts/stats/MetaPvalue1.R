library(metap)
library(rlang)
library(dplyr)
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
outdata <- read.table(
        infile, header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE
    ) %>%
    group_by(!!!syms(id_cols)) %>%
    summarise(
        N = n(),
        .pvals = list(!!sym(pval_col)),
        .groups = "drop"
    )

metaps <- c()
ns <- c()
for (ps in outdata$.pvals) {
    if (na == -1) {
        ps <- ps[!is.na(ps)]
    } else {
        ps[is.na(ps)] <- na
    }

    if (length(ps) == 0) {
        metaps <- c(metaps, NA)
        ns <- c(ns, NA)
    } else if (length(ps) == 1 && keep_single) {
        metaps <- c(metaps, ps)
        ns <- c(ns, 1)
    } else if (any(ps == 0)) {
        metaps <- c(metaps, 0)
        ns <- c(ns, length(ps))
    } else {
        metaps <- c(metaps, do.call(method, list(ps))$p)
        ns <- c(ns, length(ps))
    }
}
outdata$MetaPval <- metaps
outdata$N <- ns
outdata$.pvals <- NULL
outdata <- outdata %>% arrange(MetaPval)

if (padj != "none") {
    log$info("Calculating adjusted p-values ...")
    outdata$MetaPadj <- p.adjust(outdata$MetaPval, method = padj)

}

log$info("Writing output ...")
write.table(outdata, outfile, quote = FALSE, sep = "\t", row.names = FALSE)
