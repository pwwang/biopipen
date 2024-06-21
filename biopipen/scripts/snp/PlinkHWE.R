source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/plot.R")
library(ggprism)
theme_set(theme_prism())

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

    plotGG(
        data = hardy,
        geom = "histogram",
        outfile = paste0(output, '.hardy.png'),
        args = list(aes(x = Pval, fill = Status), alpha = 0.8, bins = 50),
        ggs = c(
            'xlab("-log10(HWE p-value)")',
            'ylab("Count")',
            'geom_vline(xintercept = -log10(cutoff), color = "red", linetype="dashed")',
            'theme(legend.position = "none")',
            'geom_text(aes(x = -log10(cutoff), y = Inf, label = cutoff), colour="red", angle=90, vjust = 1.2, hjust = 1.2)',
            'scale_fill_manual(values = c("Pass" = "blue3", "Fail" = "red3"))'  # Added line to set "Fail" color to red
        ),
        devpars = devpars
    )
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
