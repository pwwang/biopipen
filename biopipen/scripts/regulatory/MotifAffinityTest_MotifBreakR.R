library(motifbreakR)

bsgenome <- getBSgenome(genome_pkg)

# `chrom`, `start`, `end`, `name`, `score`, `strand`, `ref`, `alt`.
is_indel <- nchar(snpinfo$ref) != 1 | nchar(snpinfo$alt) != 1
snpinfo$coordname <- ifelse(
    is_indel,
    sprintf("%s:%s-%s:%s:%s", snpinfo$chrom, snpinfo$start + 1, snpinfo$end, snpinfo$ref, snpinfo$alt),
    sprintf("%s:%s:%s:%s", snpinfo$chrom, snpinfo$end, snpinfo$ref, snpinfo$alt)
)
motifbreakr_bed <- file.path(outdir, gsub("\\.vcf(\\.gz)?$|\\.bed$", ".motifbreakr.bed", basename(varfile)))
write.table(
    snpinfo[, c("chrom", "start", "end", "coordname", "score", "strand")],
    file = motifbreakr_bed,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
snps <- snps.from.file(motifbreakr_bed, search.genome = bsgenome, format = "bed", indels = any(is_indel))
snpinfo <- snpinfo[snpinfo$coordname == snps$SNP_id, , drop = FALSE]
snps@elementMetadata$SNP_id <- ifelse(
    snpinfo$name == "." | is.na(snpinfo$name) | nchar(snpinfo$name) == 0,
    snpinfo$coordname,
    snpinfo$name
)

# prepare PWMs
get_bkg <- function(base) {
    base_col <- paste0("bkg.", base)
    base_bkg <- mdb@elementMetadata[[base_col]]
    if (is.null(base_bkg) || length(base_bkg) == 0 || is.na(base_bkg[1])) {
        base_bkg <- 0.25
    } else {
        base_bkg <- as.numeric(base_bkg[1])
    }
}
bkg <- c(A = get_bkg("A"), C = get_bkg("C"), G = get_bkg("G"), T = get_bkg("T"))

# run motifbreakR
log$info("Running motifbreakR ...")
results <- motifbreakR(
    snpList = snps,
    pwmList = mdb,
    threshold = cutoff,
    method = motifbreakr_args$method,
    bkg = bkg,
    filterp = TRUE,
    show.neutral = FALSE,
    BPPARAM = MulticoreParam(ncores)
)

log$info("Calculating p values ...")
results <- calculatePvalue(results)
results$.id <- 1:length(results)
results_to_save <- as.data.frame(unname(results))
results_to_save$motifPos <- lapply(results_to_save$motifPos, function(x) paste(x, collapse = ","))
results_to_save$altPos <- lapply(results_to_save$altPos, function(x) paste(x, collapse = ","))
if (!is.null(regulator_col)) {
    results_to_save$Regulator <- in_motifs[
        match(results_to_save$providerId, in_motifs[[motif_col]]),
        regulator_col,
        drop = TRUE
    ]
}
results_to_save <- as.data.frame(apply(results_to_save, 2, as.character))

if (!is.null(motif_var_pairs)) {
    log$info("Filtering motif-variant pairs ...")
    results_to_save$motifs_vars <- paste0(results_to_save$providerId, " // ", results_to_save$SNP_id)
    results_to_save <- results_to_save[results_to_save$motifs_vars %in% motif_var_pairs, , drop = FALSE]
    results_to_save$motifs_vars <- NULL
}

write.table(
    results_to_save,
    file = file.path(outdir, "motifbreakr.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
# rm(results_to_save)

log$info("Plotting variants ...")
if (is.null(plots) || length(plots) == 0) {
    results_to_save$alleleDiff <- as.numeric(results_to_save$alleleDiff)
    results_to_save <- results_to_save[order(-abs(results_to_save$alleleDiff)), , drop = FALSE]
    results_to_save <- results_to_save[1:min(plot_nvars, nrow(results_to_save)), , drop = FALSE]
    variants <- unique(results_to_save$SNP_id)
} else {
    variants <- names(plots)
}
for (variant in variants) {
    log$info("- Variant: {variant}")
    if (is.null(plots[[variant]])) {
        plots[[variant]] <- list(devpars = devpars, which = "TRUE")
    }
    if (is.null(plots[[variant]]$which)) {
        plots[[variant]]$which <- "TRUE"
    }
    if (is.null(plots[[variant]]$devpars)) {
        plots[[variant]]$devpars <- devpars
    }
    res <- results[results$SNP_id == variant & results$.id %in% results_to_save$.id, , drop = FALSE]
    res <- subset(res, subset = eval(parse(text = plots[[variant]]$which)))

    plot_variant_motifs(res, variant, plots[[variant]]$devpars, outdir)
}
