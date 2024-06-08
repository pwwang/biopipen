source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/plot.R")
library(rlang)
library(ggprism)
theme_set(theme_prism())

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

modifier <- match.arg(modifier, c("none", "counts", "case-control", "cc", "x"))
if (modifier == "cc") { modifier <- "case-control"}

cmd <- c(
    plink,
    "--threads", ncores,
    "--bfile", input,
    "--out", output
)
if (modifier == "counts") {
    cmd <- c(cmd, "--freq", "counts")
    if (!is.list(cutoffs)) { cutoffs <- list(C1 = cutoffs) }
} else if (modifier == "case-control") {
    cmd <- c(cmd, "--freq", "case-control")
    if (!is.list(cutoffs)) { cutoffs <- list(MAF_A = cutoffs) }
} else if (modifier == "x") {
    cmd <- c(cmd, "--freqx")
    if (!is.list(cutoffs)) { cutoffs <- list("C(HOM A1)" = cutoffs) }
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

post_process <- function(suffix, snp_col = "SNP", sep = "") {
    freq <- read.table(
        paste0(output, suffix),
        header=TRUE,
        check.names=FALSE,
        row.names = NULL,
        sep = sep
    )
    iter_in <- input
    n <- 0
    for (metric_col in names(cutoffs)) {

        freq[[metric_col]] <- as.numeric(freq[[metric_col]])
        cutoff <- cutoffs[[metric_col]]
        filter <- filters[[metric_col]] %||% "no"

        if (filter == "no") {
            ge_flag <- paste0(metric_col, " >= cutoff")
            lt_flag <- paste0(metric_col, " < cutoff")
            freq$GE <- freq[[metric_col]] >= cutoff
            freq$Flag <- ifelse(freq$GE, ge_flag, lt_flag)
            freq$Flag <- factor(freq$Flag, levels = c(ge_flag, lt_flag))
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
                plotGG(
                    data = freq,
                    geom = "histogram",
                    outfile = paste0(output, suffix, ".", metric_col, ".png"),
                    args = list(aes(x = !!sym(metric_col), fill = Flag), alpha = 0.8, bins = 50),
                    ggs = c(
                        sprintf('xlab("%s")', metric_col),
                        'ylab("Count")',
                        sprintf('geom_vline(xintercept = %.3f, color = "red", linetype="dashed")', cutoff),
                        sprintf(
                            'geom_text(aes(x = %.3f, y = Inf, label = as.character(%.3f)), colour="blue", vjust = 1.5, hjust = -.1)',
                            cutoff, cutoff
                        ),
                        sprintf(
                            'scale_fill_manual(values = c("%s" = "blue3", "%s" = "green3"))',
                            ge_flag, lt_flag
                        )
                    ),
                    devpars = devpars
                )
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
            failfile <- paste0(output, suffix, ".", metric_col, ".fail")
            write.table(
                freq[[snp_col]][freq$Flag == "Fail"],
                file = failfile,
                col.names=FALSE,
                row.names=FALSE,
                quote=FALSE
            )

            if (doplot) {
                plotGG(
                    data = freq,
                    geom = "histogram",
                    outfile = paste0(output, suffix, ".", metric_col, ".png"),
                    args = list(aes(x = !!sym(metric_col), fill = Flag), alpha = 0.8, bins = 50),
                    ggs = c(
                        sprintf('xlab("%s")', metric_col),
                        'ylab("Count")',
                        sprintf('geom_vline(xintercept = %.3f, color = "blue", linetype="dashed")', cutoff),
                        sprintf(
                            'geom_text(aes(x = %.3f, y = Inf, label = as.character(%.3f)), colour="blue", vjust = 1.5, hjust = -.1)',
                            cutoff, cutoff
                        ),
                        'theme(legend.position = "none")',
                        'scale_fill_manual(values = c("Pass" = "blue3", "Fail" = "red3"))'
                    ),
                    devpars = devpars
                )
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

if (modifier == "none") {
    post_process(".frq")
} else if (modifier == "counts") {
    post_process(".frq.counts")
} else if (modifier == "case-control") {
    post_process(".frq.cc")
} else if (modifier == "x") {
    post_process(".frqx", sep = "\t")
}
