library(atSNP)
library(rtracklayer)

log_info("Converting universalmotif object to motif_library ...")

motifdb_names <- sapply(meme, function(m) m@name)
motifs <- check_motifs(motifdb_names)
meme <- filter_motifs(meme, name = motifs)
# Get the right order of motif names
motifs <- sapply(meme, function(m) m@name)

# used for atSNP
mdb <- lapply(meme, function(m) t(m@motif))
names(mdb) <- motifs

# compose one used for plotting using motifbreakR
motifdb_matrices <- lapply(meme, function(m) m@motif)
names(motifdb_matrices) <- motifs
motifdb_meta <- do.call(rbind, lapply(meme, function(m) {
    ats <- attributes(m)
    ats$dataSource <- basename(motifdb)
    ats$class <- NULL
    ats$motif <- NULL
    ats$gapinfo <- NULL
    ats$sequenceCount <- ats$nsites
    ats$providerId <- ats$name
    ats$providerName <- ats$name
    ats$organism <- if (is.null(ats$organism) || length(ats$organism) == 0) "Unknown" else ats$organism
    unlist(ats)
}))
rownames(motifdb_meta) <- motifs
pmotifs <- MotifDb:::MotifList(motifdb_matrices, tbl.metadata = motifdb_meta)

log_info("Converting snpinfo to atSNP object ...")

# c("chrom", "start", "end", "name", "score", "strand", "ref", "alt", "ref_seq", "alt_seq")
if (any(nchar(snpinfo$ref) != 1) || any(nchar(snpinfo$alt) != 1)) {
    stop("Only SNVs are supported by atSNP. Consider using motifbreakR instead if you have indels.")
}
atsnp_bed <- file.path(outdir, gsub("\\.vcf(\\.gz)?$|\\.bed$", ".atsnp.txt", basename(varfile)))
snpinfo$name <- ifelse(
    snpinfo$name == "." | is.na(snpinfo$name) | nchar(snpinfo$name) == 0,
    sprintf("%s:%s", snpinfo$chrom, snpinfo$end),
    snpinfo$name
)
snpinfo$a1 <- snpinfo$ref
snpinfo$a2 <- snpinfo$alt
snpinfo$chr <- snpinfo$chrom
snpinfo$snp <- snpinfo$end
snpinfo$snpid <- snpinfo$name
write.table(
    snpinfo[, c("snpid", "a1", "a2", "chr", "snp")],
    file = atsnp_bed,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)
k <- max(sapply(mdb, nrow))
snps <- LoadSNPData(
    atsnp_bed,
    genome.lib = genome_pkg,
    mutation = TRUE,  # force using given ref and alt
    default.par = nrow(snpinfo) < 1000,
    half.window.size = k
)

# run motifbreakR
log_info("Running atSNP ...")
atsnp_scores <- ComputeMotifScore(mdb, snps, ncores = ncores)

log_info("Calculating p values ...")
atsnp_result <- ComputePValues(
    motif.lib = mdb,
    snp.info = snps,
    motif.scores = atsnp_scores$motif.scores,
    ncores = ncores,
    testing.mc = TRUE
)

padj_col <- paste0(atsnp_args$p, "_adj")
atsnp_result[[padj_col]] <- p.adjust(atsnp_result[[atsnp_args$p]], method = atsnp_args$padj)
cutoff_col <- if (atsnp_args$padj_cutoff) padj_col else atsnp_args$p
atsnp_result <- atsnp_result[atsnp_result[[cutoff_col]] < cutoff, , drop = FALSE]
# order by p value
atsnp_result <- atsnp_result[order(atsnp_result[[cutoff_col]]), , drop = FALSE]
snpinfo <- snpinfo[match(atsnp_result$snpid, snpinfo$snpid), , drop = FALSE]
atsnp_result$chr <- snpinfo$chr
atsnp_result$start <- snpinfo$start
atsnp_result$end <- snpinfo$end
atsnp_result$SNP_id <- snpinfo$snpid
atsnp_result$snpid <- NULL
atsnp_result$REF <- snpinfo$ref
atsnp_result$ALT <- snpinfo$alt
atsnp_result$providerName <- atsnp_result$motif
atsnp_result$providerId <- atsnp_result$providerName <- atsnp_result$motif
atsnp_result$motif <- NULL
atsnp_result$strand <- snpinfo$strand
atsnp_result$score <- snpinfo$score
atsnp_result$snpbase <- NULL
atsnp_result$altPos <- 1
atsnp_result$varType <- "SNV"
atsnp_result$motifPos <- sapply(1:nrow(atsnp_result), function(i) {
    paste(c(atsnp_result$ref_start[i] - k, atsnp_result$ref_end[i] - k), collapse = ",")
})
if (!is.null(regulator_col)) {
    atsnp_result$Regulator <- in_motifs[
        match(atsnp_result$providerId, in_motifs[[motif_col]]),
        regulator_col,
        drop = TRUE
    ]
}

write.table(
    atsnp_result,
    file = file.path(outdir, "atsnp.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

log_info("Plotting variants ...")
# Convert result to GRanges object
atsnp_result$alleleDiff <- -atsnp_result[[cutoff_col]]
atsnp_result$effect <- "strong"
atsnp_result$motifPos <- lapply(atsnp_result$motifPos, function(x) as.integer(unlist(strsplit(x, ","))))
atsnp_result <- makeGRangesFromDataFrame(atsnp_result, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
attributes(atsnp_result)$genome.package <- genome_pkg
attributes(atsnp_result)$motifs <- pmotifs

plot_variant(atsnp_result)
