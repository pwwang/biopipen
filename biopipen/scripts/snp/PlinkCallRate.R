library(plotthis)
library(biopipen.utils)

indir <- {{in.indir | r}}
outdir <- {{out.outdir | r}}
plink <- {{envs.plink | r}}
ncores <- {{envs.ncores | r}}
doplot <- {{envs.plot | r}}
devpars <- {{envs.devpars | r}}
samplecr <- {{envs.samplecr | r}}
varcr <- {{envs.varcr | r}}
max_iter <- {{envs.max_iter | r}}

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

all_smiss_file = paste0(output, '.smiss')
all_vmiss_file = paste0(output, '.vmiss')
all_samplecr_fail_file = paste0(output, '.samplecr.fail')
all_varcr_fail_file = paste0(output, '.varcr.fail')
if (file.exists(all_smiss_file)) invisible(file.remove(all_smiss_file))
if (file.exists(all_vmiss_file)) invisible(file.remove(all_vmiss_file))
for (i in 1:max_iter) {
    log$info("Iteration {i} ...")
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

    smissfile <- paste0(iter_out, '.smiss')
    smiss <- read.table(
        smissfile,
        header = TRUE,
        row.names = NULL,
        check.names = FALSE,
        comment.char = ""
    )
    smiss$Iteration <- i
    # append it to all_smiss_file
    write.table(
        smiss,
        all_smiss_file,
        append = i > 1,
        col.names = !file.exists(all_smiss_file),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
    callrate.sample <- data.frame(Callrate = 1 - smiss$F_MISS)
    rownames(callrate.sample) <- paste(smiss$FID, smiss$IID, sep = "\t")
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

    vmiss <- read.table(
        paste0(iter_out, '.vmiss'),
        header = TRUE,
        row.names = NULL,
        check.names = FALSE,
        comment.char = ""
    )
    vmiss$Iteration <- i
    # append it to all_vmiss_file
    write.table(
        vmiss,
        all_vmiss_file,
        append = i > 1,
        col.names = !file.exists(all_vmiss_file),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
    vmiss$Callrate <- 1 - vmiss$F_MISS
    callrate.var.fail <- vmiss[which(vmiss$Callrate < varcr), 'ID', drop = TRUE]
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

smiss <- read.table(
    smissfile,
    header = TRUE,
    row.names = NULL,
    check.names = FALSE,
    comment.char = ""
)
callrate.sample <- data.frame(Callrate = 1 - smiss$F_MISS)
rownames(callrate.sample) <- paste(smiss$FID, smiss$IID, sep = "\t")

vmiss <- read.table(
    paste0(iter_out, '.vmiss'),
    header = TRUE,
    row.names = NULL,
    check.names = FALSE,
    comment.char = ""
)
vmiss$Callrate <- 1 - vmiss$F_MISS

if (doplot) {
    log$info("Plotting ...")
    callrate.sample$Status <- "Pass"
    callrate.sample[callrate.sample.fail, "Status"] <- "Fail"
    callrate.sample$Status <- factor(callrate.sample$Status, levels = c("Fail", "Pass"))

    p_callrate_file <- paste0(output, '.samplecr.png')
    p_callrate <- Histogram(
        callrate.sample,
        x = "Callrate",
        group_by = "Status",
        xlab = "Sample Call Rate",
        ylab = "Count",
        palette = "Set1",
        alpha = 0.8,
        bins = 50
    )
    res <- 70
    height <- attr(p_callrate, "height") * res
    width <- attr(p_callrate, "width") * res
    png(p_callrate_file, width = width, height = height, res = res)
    print(p_callrate)
    dev.off()

    vmiss$Status <- "Pass"
    vmiss[which(vmiss$Callrate < varcr), "Status"] <- "Fail"
    vmiss$Status <- factor(vmiss$Status, levels = c("Fail", "Pass"))

    p_varcr_file <- paste0(output, '.varcr.png')
    p_varcr <- Histogram(
        vmiss,
        x = "Callrate",
        group_by = "Status",
        xlab = "Variant Call Rate",
        ylab = "Count",
        palette = "Set1",
        alpha = 0.8,
        bins = 50
    )
    res <- 70
    height <- attr(p_varcr, "height") * res
    width <- attr(p_varcr, "width") * res
    png(p_varcr_file, width = width, height = height, res = res)
    print(p_varcr)
    dev.off()
}
