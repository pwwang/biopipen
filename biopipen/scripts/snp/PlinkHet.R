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

    plotGG(
        data = het,
        geom = "histogram",
        outfile = paste0(output, '.het.png'),
        args = list(aes(fill = Status, x = Het), alpha = 0.8, bins = 50),
        ggs = c(
            'xlab("Sample Heterozygosity")',
            'ylab("Count")',
            'geom_vline(xintercept = c(het.mean-cutoff*het.sd, het.mean+cutoff*het.sd), color = "red", linetype="dashed")',
            'geom_vline(xintercept = het.mean, color = "blue", linetype="dashed")',
            'theme(legend.position = "none")',
            'geom_text(aes(x = het.mean-cutoff*het.sd, y = Inf, label = sprintf("mean - %ssd (%.3f)", cutoff, het.mean - cutoff*het.sd)), colour="red", angle=90, vjust = 1.2, hjust = 1.2)',
            'geom_text(aes(x = het.mean+cutoff*het.sd, y = Inf, label = sprintf("mean + %ssd (%.3f)", cutoff, het.mean + cutoff*het.sd)), colour="red", angle=90, vjust = 1.2, hjust = 1.2)',
            'geom_text(aes(x = het.mean, y = Inf, label = sprintf("mean (%.3f)", het.mean)), colour="blue", vjust = 1.5, hjust = -.1)',
            'scale_fill_manual(values = c("Pass" = "blue3", "Fail" = "red3"))'
        ),
        devpars = devpars
    )
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
