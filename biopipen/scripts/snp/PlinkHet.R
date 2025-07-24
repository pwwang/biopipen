library(plotthis)
library(biopipen.utils)

indir <- {{in.indir | r}}
outdir <- {{out.outdir | r}}
plink <- {{envs.plink | r}}
ncores <- {{envs.ncores | r}}
cutoff <- {{envs.cutoff | r}}
doplot <- {{envs.plot | r}}
devpars <- {{envs.devpars | r}}

log <- get_logger()

bedfile = Sys.glob(file.path(indir, '*.bed'))
if (length(bedfile) == 0)
    stop("No bed files found in the input directory.")
if (length(bedfile) > 1) {
    log$warn("Multiple bed files found in the input directory. Using the first one.")
    bedfile <- bedfile[1]
}
input <- tools::file_path_sans_ext(bedfile)
output <- file.path(outdir, basename(input))

# need .afreq for --het for plink2
freq_cmd <- cmd <- c(
    plink,
    "--threads", ncores,
    "--bfile", input,
    "--freq",
    "--out", output
)
run_command(freq_cmd, fg = TRUE)

cmd <- c(
    plink,
    "--threads", ncores,
    "--bfile", input,
    "--het",
    "--out", output,
    "--read-freq", paste0(output, '.afreq')
)
run_command(cmd, fg = TRUE)

phet <- read.table(
    paste0(output, '.het'),
    header = TRUE,
    row.names = NULL,
    check.names = FALSE,
    comment.char = ""
)
het <- data.frame(Het = 1 - phet[, "O(HOM)"]/phet[, "OBS_CT"])
rownames(het) <- paste(phet$FID, phet$IID, sep = "\t")
het.mean <- mean(het$Het, na.rm = TRUE)
het.sd <- sd(het$Het, na.rm = TRUE)
het.fail <- rownames(het[
    !is.na(het$Het) & (het$Het < het.mean-cutoff*het.sd | het$Het > het.mean+cutoff*het.sd), , drop = FALSE
])
writeLines(het.fail, con = file(paste0(output, '.het.fail')))

if (doplot) {
    het$Status <- "Pass"
    het[het.fail, "Status"] <- "Fail"
    het$Status <- factor(het$Status, levels = c("Fail", "Pass"))

    p <- Histogram(
        het,
        x = "Het",
        group_by = "Status",
        alpha = 0.8,
        bins = 50,
        xlab = "Sample Heterozygosity",
        ylab = "Count",
        palette = "Set1"
    )
    res <- 70
    height <- attr(p, "height") * res
    width <- attr(p, "width") * res
    png(
        filename = paste0(output, '.het.png'),
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
    "--remove", paste0(output, '.het.fail'),
    "--make-bed",
    "--out", output
)
run_command(cmd, fg = TRUE)
