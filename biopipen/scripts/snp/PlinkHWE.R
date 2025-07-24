library(plotthis)
library(biopipen.utils)

indir <- {{in.indir | r}}
outdir <- {{out.outdir | r}}
plink <- {{envs.plink | r}}
ncores <- {{envs.ncores | r}}
cutoff <- {{envs.cutoff | r}}
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

cmd <- c(
    plink,
    "--threads", ncores,
    "--bfile", input,
    "--hardy",
    "--out", output
)
run_command(cmd, fg = TRUE)

hardy <- read.table(
    paste0(output, '.hardy'),
    header = TRUE,
    row.names = NULL,
    check.names = FALSE,
    comment.char = ""
)
hardy.fail <- hardy[which(hardy$P < cutoff), 'ID', drop = FALSE]
write.table(
    hardy.fail,
    paste0(output, '.hardy.fail'),
    col.names = FALSE,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
)

if (doplot) {
    hardy$Pval <- -log10(hardy$P)
    hardy$Status <- "Pass"
    hardy[which(hardy$SNP %in% hardy.fail$SNP), "Status"] <- "Fail"
    hardy$Status <- factor(hardy$Status, levels = c("Fail", "Pass"))

    p <- Histogram(
        hardy,
        x = "Pval",
        group_by = "Status",
        alpha = 0.8,
        bins = 50,
        xlab = "-log10(HWE p-value)",
        ylab = "Count",
        palette = "Set1"
    )
    res <- 70
    height <- attr(p, "height") * res
    width <- attr(p, "width") * res
    png(
        filename = paste0(output, '.hardy.png'),
        width = width,
        height = height,
        res = res
    )
    print(p)
    dev.off()
}

cmd <- c(
    plink,
    "--threads", ncores,
    "--bfile", input,
    "--exclude", paste0(output, '.hardy.fail'),
    "--make-bed",
    "--out", output
)
run_command(cmd, fg = TRUE)
