library(motifbreakR)
bsgenome <- getBSgenome(genome_pkg)

log_info("Converting universalmotif object to MotifDb object ...")

motifdb_names <- sapply(meme, function(m) m@name)
motifs <- check_motifs(motifdb_names)
meme <- filter_motifs(meme, name = motifs)
# Get the right order of motif names
motifs <- sapply(meme, function(m) m@name)
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
mdb <- MotifDb:::MotifList(motifdb_matrices, tbl.metadata = motifdb_meta)

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
log_info("Running motifbreakR ...")
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

log_info("Calculating p values ...")
results <- calculatePvalue(results)
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
results_to_save <- apply(results_to_save, 2, as.character)

write.table(
    results_to_save,
    file = file.path(outdir, "motifbreakr.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
rm(results_to_save)

plot_variant(results)
