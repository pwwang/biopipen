library(metap)
library(rlang)
library(dplyr)
library(biopipen.utils)

infiles <- {{in.infiles | each: str | r}}
outfile <- {{out.outfile | r}}
id_cols <- {{envs.id_cols | r}}
id_exprs <- {{envs.id_exprs | r}}
pval_cols <- {{envs.pval_cols | r}}
method <- {{envs.method | r}}
na <- {{envs.na | r}}
keep_single <- {{envs.keep_single | r}}
padj <- {{envs.padj | r}}

if (method == "fisher") { method = "sumlog" }

log <- get_logger()

if (length(infiles) == 1 && padj == "none") {
    log$info("Only one input file, copying to output ...")
    file.copy(infiles, outfile)
} else if (length(infiles) == 1) {
    log$info("Only one input file, performing p-value adjustment ...")
    if (is.null(pval_cols)) {
        stop("Must provide envs.pval_cols")
    }
    indata <- read.table(infiles, header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE)
    if (!pval_cols %in% colnames(indata)) {
        stop("envs.pval_cols does not exist in input file")
    }
    indata$Padj <- p.adjust(indata[, pval_cols], method = padj)

    log$info("Writing output ...")
    write.table(indata, outfile, quote = FALSE, sep = "\t", row.names = FALSE)
} else {
    # Check pval_cols
    if (is.null(pval_cols)) {
        stop("Must provide envs.pval_cols")
    }
    if (length(pval_cols) == 1) {
        pval_cols <- trimws(strsplit(pval_cols, ",")[[1]])
    }
    if (length(pval_cols) == 1) {
        pval_cols <- rep(pval_cols, length(infiles))
    }
    if (length(pval_cols) != length(infiles)) {
        stop("envs.pval_cols must be a single name or have the same length as in.infiles")
    }

    # Check id_cols
    if (is.null(id_cols)) {
        stop("Must provide envs.id_cols")
    }
    if (length(id_cols) == 1) {
        id_cols <- trimws(strsplit(id_cols, ",")[[1]])
    }

    # Check id_exprs
    if (!is.null(id_exprs)) {
        if (length(id_exprs) == 1) {
            id_exprs <- rep(id_exprs, length(infiles))
        }
        if (length(id_exprs) != length(infiles)) {
            stop("envs.id_exprs must be a single expression or have the same length as in.infiles")
        }
        if (length(id_cols) != 1) {
            stop("envs.id_cols must be a single name if envs.id_exprs is provided")
        }
    }

    log$info("Reading and preparing data ...")
    outdata <- NULL
    for (i in seq_along(infiles)) {
        infile <- infiles[i]
        name <- tools::file_path_sans_ext(basename(infile))
        pval_col <- paste0("Pval_", name)
        dat <- read.table(
            infile, header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE
        )
        if (!is.null(id_exprs)) {
            dat <- dat %>% mutate(!!sym(id_cols) := !!parse_expr(id_exprs[i]))
        }
        dat <- dat %>% dplyr::select(all_of(id_cols), !!sym(pval_col) := !!sym(pval_cols[i]))

        if (is.null(outdata)) {
            outdata <- dat
        } else {
            outdata <- full_join(outdata, dat, by = id_cols)
        }
    }

    log$info("Running metap on each row ...")
    metaps <- c()
    ns <- c()
    pval_columns <- setdiff(colnames(outdata), id_cols)
    for (i in seq_len(nrow(outdata))) {
        ps <- unlist(outdata[i, pval_columns, drop = TRUE])
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
    outdata <- outdata %>% arrange(MetaPval)

    if (padj != "none") {
        log$info("Calculating adjusted p-values ...")
        outdata$MetaPadj <- p.adjust(outdata$MetaPval, method = padj)

    }

    log$info("Writing output ...")
    write.table(outdata, outfile, quote = FALSE, sep = "\t", row.names = FALSE)
}
