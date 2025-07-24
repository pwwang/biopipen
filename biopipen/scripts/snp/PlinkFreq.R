library(rlang)
library(plotthis)
library(biopipen.utils)

indir <- {{in.indir | r}}
outdir <- {{out.outdir | r}}
plink <- {{envs.plink | r}}
ncores <- {{envs.ncores | r}}
modifier <- {{envs.modifier | r}}
gz <- {{envs.gz | r}}
cutoffs <- {{envs.cutoff | r}}
filters <- {{envs.filter | r}}
doplot <- {{envs.plot | r}}
devpars <- {{envs.devpars | r}}

bedfile = Sys.glob(file.path(indir, '*.bed'))
if (length(bedfile) == 0)
    stop("No bed files found in the input directory.")
if (length(bedfile) > 1) {
    log_warn("Multiple bed files found in the input directory. Using the first one.")
    bedfile <- bedfile[1]
}
input <- tools::file_path_sans_ext(bedfile)
output <- file.path(outdir, basename(input))

modifier <- match.arg(modifier, c("none", "counts", "x"))

cmd <- c(
    plink,
    "--threads", ncores,
    "--bfile", input,
    "--out", output
)
if (modifier == "counts") {
    cmd <- c(cmd, "--freq", "counts")
    if (!is.list(cutoffs)) { cutoffs <- list(ALT1_CT = cutoffs) }
# } else if (modifier == "case-control") {
#     cmd <- c(cmd, "--freq", "case-control")
#     if (!is.list(cutoffs)) { cutoffs <- list(MAF_A = cutoffs) }
} else if (modifier == "x") {
    cmd <- c(cmd, "--geno-counts")
    if (!is.list(cutoffs)) { cutoffs <- list("HOM_ALT1_CT" = cutoffs) }
} else {
    cmd <- c(cmd, "--freq")
    if (!is.list(cutoffs)) { cutoffs <- list(MAF = cutoffs) }
}
if (isTRUE(gz)) { cmd <- c(cmd, "gz") }

if (!is.list(filters)) {
    filters <- as.list(rep(filters, length(cutoffs)))
    names(filters) <- names(cutoffs)
} else {
    for (name in names(filters)) {
        if (is.null(cutoffs[[name]])) {
            stop(paste0("Cutoff for filter ", name, " is not provided."))
        }
    }
}

run_command(cmd, fg = TRUE)

post_process <- function(suffix, snp_col = "ID", sep = "\t", modifier = NULL) {
    freq <- read.table(
        paste0(output, suffix),
        header=TRUE,
        check.names=FALSE,
        row.names = NULL,
        sep = sep,
        comment = ""
    )
    colnames(freq)[1] <- sub("#", "", colnames(freq)[1])
    if (!is.null(modifier)) { freq <- modifier(freq) }
    iter_in <- input
    n <- 0
    for (metric_col in names(cutoffs)) {
        if (is.null(cutoffs[[metric_col]])) {
            stop(paste0(
                "Cutoff for metric ",
                metric_col,
                " is not provided in ",
                suffix, "(x) file."))
        }

        freq[[metric_col]] <- as.numeric(freq[[metric_col]])
        cutoff <- cutoffs[[metric_col]]
        filter <- filters[[metric_col]] %||% "no"

        if (filter == "no") {
            ge_flag <- paste0(metric_col, " >= ", cutoff)
            lt_flag <- paste0(metric_col, " < ", cutoff)
            freq$GE <- freq[[metric_col]] >= cutoff
            freq$Flag <- ifelse(freq$GE, ge_flag, lt_flag)
            freq$Flag <- factor(freq$Flag, levels = c(lt_flag, ge_flag))
            write.table(
                freq[[snp_col]][freq$GE],
                file = paste0(output, suffix, ".", metric_col, ".ge"),
                col.names=FALSE,
                row.names=FALSE,
                quote=FALSE
            )
            write.table(
                freq[[snp_col]][!freq$GE],
                file = paste0(output, suffix, ".", metric_col, ".lt"),
                col.names=FALSE,
                row.names=FALSE,
                quote=FALSE
            )

            if (doplot) {
                p <- Histogram(
                    freq,
                    x = metric_col,
                    group_by = "Flag",
                    alpha = 0.8,
                    bins = 50,
                    xlab = metric_col,
                    ylab = "Count",
                    palette = "Set1"
                )
                res <- 70
                height <- attr(p, "height") * res
                width <- attr(p, "width") * res
                png(paste0(output, suffix, ".", metric_col, ".png"), width = width, height = height, res = res)
                print(p)
                dev.off()
            }
        } else {
            iter_dir <- file.path(outdir, paste0(metric_col, "_filtered"))
            dir.create(iter_dir, showWarnings = FALSE)
            iter_out <- file.path(iter_dir, basename(output))

            filter <- match.arg(filter, c("gt", "lt", "ge", "le"))
            indicate <- function(metric){
                if (filter == "gt") {
                    return(freq[[metric_col]] > cutoff)
                } else if (filter == "lt") {
                    return(freq[[metric_col]] < cutoff)
                } else if (filter == "ge") {
                    return(freq[[metric_col]] >= cutoff)
                } else if (filter == "le") {
                    return(freq[[metric_col]] <= cutoff)
                }
            }
            freq$Flag <- ifelse(indicate(freq), "Fail", "Pass")
            freq$Flag <- factor(freq$Flag, levels = c("Fail", "Pass"))
            failfile <- paste0(output, suffix, ".", metric_col, ".fail")
            write.table(
                freq[[snp_col]][freq$Flag == "Fail"],
                file = failfile,
                col.names=FALSE,
                row.names=FALSE,
                quote=FALSE
            )

            if (doplot) {
                p <- Histogram(
                    freq,
                    x = metric_col,
                    group_by = "Flag",
                    alpha = 0.8,
                    bins = 50,
                    xlab = metric_col,
                    ylab = "Count",
                    palette = "Set1"
                )
                res <- 70
                height <- attr(p, "height") * res
                width <- attr(p, "width") * res
                png(paste0(output, suffix, ".", metric_col, ".png"), width = width, height = height, res = res)
                print(p)
                dev.off()
            }

            filter_cmd <- c(
                plink,
                "--threads", ncores,
                "--bfile", shQuote(iter_in),
                "--exclude", shQuote(failfile),
                "--make-bed",
                "--out", shQuote(iter_out)
            )
            run_command(filter_cmd, fg = TRUE)

            iter_in <- iter_out
            n <- n + 1

            if (n == length(cutoffs)) {
                # make symbolic links to output from input .bed, .bim and .fam files
                file.symlink(paste0(iter_in, '.bed'), paste0(output, '.bed'))
                file.symlink(paste0(iter_in, '.bim'), paste0(output, '.bim'))
                file.symlink(paste0(iter_in, '.fam'), paste0(output, '.fam'))
            }
        }
    }
}

splitup <- function(x, agg = NULL) {
    sp <- strsplit(as.character(x), ",")
    if (is.null(agg)) {
        return(sp)
    }
    return(sapply(sp, agg))
}
if (modifier == "none") {
    mod <- function(freq) {
        # Add ALT1, ALT1_FREQ, REF_FREQ and MAF columns
        writing = FALSE
        if (is.null(freq$ALT1)) {
            # should be the first allele of ALT
            freq$ALT1 <- splitup(freq$ALT, agg = function(s) s[1])
            writing = TRUE
        }
        if (is.null(freq$ALT1_FREQ)) {
            freq$ALT1_FREQ <- as.double(splitup(freq$ALT_FREQS, agg = function(s) s[1]))
            writing = TRUE
        }
        if (is.null(freq$REF_FREQ)) {
            freq$REF_FREQ <- 1 - splitup(freq$ALT_FREQS, agg = function(s) sum(as.double(s)))
            writing = TRUE
        }
        if (is.null(freq$MAF)) {
            min_alt_freqs <- splitup(freq$ALT_FREQS, agg = function(s) min(as.double(s)))
            freq$MAF <- pmin(freq$REF_FREQ, min_alt_freqs)
            writing = TRUE
        }
        if (writing) {
            write.table(
                freq,
                file = paste0(output, ".afreqx"),
                col.names=TRUE,
                row.names=FALSE,
                quote=FALSE,
                sep = "\t"
            )
        }
        return(freq)
    }
    post_process(".afreq", modifier = mod)
} else if (modifier == "counts") {
    mod <- function(freq) {
        # Add ALT1, ALT1_CT, and REF_CT columns
        writing = FALSE
        if (is.null(freq$ALT1)) {
            # should be the first allele of ALT
            freq$ALT1 <- splitup(freq$ALT, agg = function(s) s[1])
            writing = TRUE
        }
        if (is.null(freq$ALT1_CT)) {
            freq$ALT1_CT <- as.integer(splitup(freq$ALT_CTS, agg = function(s) s[1]))
            writing = TRUE
        }
        if (is.null(freq$REF_CT)) {
            freq$REF_CT <- freq$OBS_CT - splitup(freq$ALT_CTS, agg = function(s) sum(as.integer(s)))
            writing = TRUE
        }
        if (writing) {
            write.table(
                freq,
                file = paste0(output, ".acountx"),
                col.names=TRUE,
                row.names=FALSE,
                quote=FALSE,
                sep = "\t"
            )
        }
        return(freq)
    }
    post_process(".acount", modifier = mod)
# } else if (modifier == "case-control") {
#     post_process(".frq.cc")
} else if (modifier == "x") {
    mod <- function(freq) {
        # Add ALT1, HET_REF_ALT1_CT, HOM_ALT1_CT
        writing = FALSE
        if (is.null(freq$ALT1)) {
            # should be the first allele of ALT
            freq$ALT1 <- splitup(freq$ALT, agg = function(s) s[1])
            writing = TRUE
        }
        if (is.null(freq$HET_REF_ALT1_CT)) {
            freq$HET_REF_ALT1_CT <- as.integer(splitup(freq$HET_REF_ALT_CTS, agg = function(s) s[1]))
            writing = TRUE
        }
        if (is.null(freq$HOM_ALT1_CT)) {
            freq$HOM_ALT1_CT <- as.integer(splitup(freq$TWO_ALT_GENO_CTS, agg = function(s) s[1]))
            writing = TRUE
        }
        return(freq)
    }
    post_process(".gcount", modifier = mod)
}
