source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/plot.R")
library(ggprism)
theme_set(theme_prism())

indir <- {{in.indir | r}}
outdir <- {{out.outdir | r}}
plink <- {{envs.plink | r}}
ncores <- {{envs.ncores | r}}
doplot <- {{envs.plot | r}}
devpars <- {{envs.devpars | r}}
samplecr <- {{envs.samplecr | r}}
varcr <- {{envs.varcr | r}}
max_iter <- {{envs.max_iter | r}}

bedfile = Sys.glob(file.path(indir, '*.bed'))
if (length(bedfile) == 0)
    stop("No bed files found in the input directory.")
if (length(bedfile) > 1) {
    log_warn("Multiple bed files found in the input directory. Using the first one.")
    bedfile <- bedfile[1]
}
input <- tools::file_path_sans_ext(bedfile)
output <- file.path(outdir, basename(input))

all_imiss_file = paste0(output, '.imiss')
all_lmiss_file = paste0(output, '.lmiss')
all_samplecr_fail_file = paste0(output, '.samplecr.fail')
all_varcr_fail_file = paste0(output, '.varcr.fail')
for (i in 1:max_iter) {
    log_info("Iteration {i} ...")
    # iter_out <- paste0(output, "-", i)
    iter_dir <- file.path(outdir, paste0("iter", i))
    dir.create(iter_dir, showWarnings = FALSE)
    iter_out <- file.path(iter_dir, basename(output))
    cmd <- c(
        plink,
        "--threads", ncores,
        "--bfile", input,
        "--missing",
        "--out", iter_out
    )
    run_command(cmd, fg = TRUE)

    imissfile <- paste0(iter_out, '.imiss')
    imiss <- read.table(
        imissfile,
        header = TRUE,
        row.names = NULL,
        check.names = FALSE
    )
    imiss$Iteration <- i
    # append it to all_imiss_file
    write.table(
        imiss,
        all_imiss_file,
        append = i > 1,
        col.names = !file.exists(all_imiss_file),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
    callrate.sample <- data.frame(Callrate = 1 - imiss$F_MISS)
    rownames(callrate.sample) <- paste(imiss$FID, imiss$IID, sep = "\t")
    callrate.sample.fail = rownames(callrate.sample[
        callrate.sample$Callrate < samplecr, , drop = FALSE
    ])
    writeLines(callrate.sample.fail, con = file(paste0(iter_out, '.samplecr.fail')))
    # append it to all_samplecr_fail_file
    write(
        paste0(sapply(
            callrate.sample.fail,
            function(x){ paste0(x, "\n") }
        ), collapse = ""),
        file = file(all_samplecr_fail_file),
        append = i > 1
    )

    lmiss <- read.table(
        paste0(iter_out, '.lmiss'),
        header = TRUE,
        row.names = NULL,
        check.names = FALSE
    )
    lmiss$Iteration <- i
    # append it to all_lmiss_file
    write.table(
        lmiss,
        all_lmiss_file,
        append = i > 1,
        col.names = !file.exists(all_lmiss_file),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
    lmiss$Callrate <- 1 - lmiss$F_MISS
    callrate.var.fail <- lmiss[which(lmiss$Callrate < varcr), 'SNP', drop = TRUE]
    writeLines(callrate.var.fail, con = file(paste0(iter_out, '.varcr.fail')))
    # append it to all_varcr_fail_file
    write(
        paste0(sapply(
            callrate.var.fail,
            function(x){ paste0(x, "\n") }
        ), collapse = ""),
        file = file(all_varcr_fail_file),
        append = i > 1
    )

    if (length(callrate.sample.fail) == 0 && length(callrate.var.fail) == 0) {
        # make symbolic links to output from input .bed, .bim and .fam files
        file.symlink(paste0(input, '.bed'), paste0(output, '.bed'))
        file.symlink(paste0(input, '.bim'), paste0(output, '.bim'))
        file.symlink(paste0(input, '.fam'), paste0(output, '.fam'))
        break
    }

    # remove samples in iter_out.samplecr.fail and variants in iter_out.varcr.fail
    cmd <- c(
        plink,
        "--threads", ncores,
        "--bfile", input,
        "--remove", paste0(iter_out, '.samplecr.fail'),
        "--exclude", paste0(iter_out, '.varcr.fail'),
        "--make-bed",
        "--out", iter_out
    )
    run_command(cmd, fg = TRUE)
    input <- iter_out
}

imiss <- read.table(
    imissfile,
    header = TRUE,
    row.names = NULL,
    check.names = FALSE
)
callrate.sample <- data.frame(Callrate = 1 - imiss$F_MISS)
rownames(callrate.sample) <- paste(imiss$FID, imiss$IID, sep = "\t")

lmiss <- read.table(
    paste0(iter_out, '.lmiss'),
    header = TRUE,
    row.names = NULL,
    check.names = FALSE
)
lmiss$Callrate <- 1 - lmiss$F_MISS

if (doplot) {
    log_info("Plotting ...")
    callrate.sample$Status <- "Pass"
    callrate.sample[callrate.sample.fail, "Status"] <- "Fail"
    plotGG(
        data = callrate.sample,
        geom = "histogram",
        outfile = paste0(output, '.samplecr.png'),
        args = list(aes(fill = Status, x = Callrate), alpha = 0.8, bins = 50),
        ggs = c(
            'xlab("Sample Call Rate")',
            'ylab("Count")',
            'geom_vline(xintercept = samplecr, color = "red", linetype="dashed")',
            'theme(legend.position = "none")',
            'geom_text(aes(x = samplecr, y = Inf, label = samplecr), colour="red", angle=90, vjust = 1.2, hjust = 1.2)',
            'scale_fill_manual(values = c("Pass" = "blue3", "Fail" = "red3"))'
        )
    )

    lmiss$Status <- "Pass"
    lmiss[which(lmiss$Callrate < varcr), "Status"] <- "Fail"
    plotGG(
        data = lmiss,
        geom = "histogram",
        outfile = paste0(output, '.varcr.png'),
        args = list(aes(fill = Status, x = Callrate), alpha = 0.8, bins = 50),
        ggs = c(
            'xlab("Variant Call Rate")',
            'ylab("Count")',
            'geom_vline(xintercept = varcr, color = "red", linetype="dashed")',
            'theme(legend.position = "none")',
            'geom_text(aes(x = varcr, y = Inf, label = varcr), colour="red", angle=90, vjust = 1.2, hjust = 1.2)',
            'scale_fill_manual(values = c("Pass" = "blue3", "Fail" = "red3"))'
        ),
        devpars = devpars
    )
}
