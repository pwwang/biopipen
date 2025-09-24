{% include biopipen_dir + "/scripts/regulatory/motifs-common.R" %}

library(BSgenome)
library(GenomicRanges)
library(biopipen.utils)

infile <- {{in.infile | r}}
outdir <- {{out.outdir | r}}
genome <- {{envs.genome | r}}
motifdb <- {{envs.motifdb | r}}
motif_col <- {{envs.motif_col | r}}
regulator_col <- {{envs.regulator_col | r}}
regmotifs <- {{envs.regmotifs | r}}
notfound <- {{envs.notfound | r}}
devpars <- {{envs.devpars | r}}
plot_vars <- {{envs.plot_vars | r}}

if (is.null(motifdb) || !file.exists(motifdb)) {
    stop("Motif database (envs.motifdb) is required and must exist")
}

if (is.null(genome)) {
    stop("Reference genome (envs.ref) is required and must exist")
}

if (is.null(motif_col) && is.null(regulator_col)) {
    stop("Either motif (envs.motif_col) or regulator (envs.regulator_col) column must be provided")
}

log <- get_logger()

log$info("Reading input data ...")
indata <- read.table(infile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)

log$info("Ensuring regulators in the input data ...")
indata <- ensure_regulator_motifs(indata, outdir, motif_col, regulator_col, "SNP_id", regmotifs, notfound = notfound)
genome_pkg <- get_genome_pkg(genome)

log$info("Reading motif database ...")
meme <- read_meme_to_motifdb(motifdb, indata, motif_col, regulator_col, notfound, outdir)

log$info("Composing motifbreakR results from input data ...")
indata$chr <- indata$chrom %||% indata$chr %||% indata$seqnames
indata$seqnames <- NULL
indata$strand <- indata$strand %||% "+"
indata$varType <- indata$varType %||% "SNV"
indata$geneSymbol <- indata$geneSymbol %||% indata$Regulator
indata$providerId <- indata$providerId %||% indata$motif
indata$providerName <- indata$providerName %||% indata$providerId
indata$dataSource <- indata$dataSource %||% strsplit(basename(motifdb), "\\.")[[1]][1]
indata$effect <- indata$effect %||% "strong"
indata$altPos <- indata$altPos %||% 1
indata$alleleDiff <- indata$alleleDiff %||% indata$score %||% 0

# check other required columns
for (col in c("start", "end", "SNP_id", "REF", "ALT", "motifPos")) {
    if (!(col %in% colnames(indata))) {
        stop("Column '", col, "' is required in the input data")
    }
}
indata$motifPos <- lapply(indata$motifPos, function(x) as.integer(unlist(strsplit(x, ","))))
indata <- makeGRangesFromDataFrame(indata, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
genome(indata) <- genome
attributes(indata)$genome.package <- genome_pkg
attributes(indata)$motifs <- meme

log$info("Plotting variants ...")
if (is.null(plot_vars)) {
    plot_vars <- unique(indata$SNP_id)
} else if (length(plot_vars) > 1) {
    plot_vars <- unique(plot_vars)
} else {
    plot_vars <- strsplit(plot_vars, ",")[[1]]
}
for (pvar in plot_vars) {
    log$info("- Variant: {pvar}")
    plot_variant_motifs(indata, pvar, devpars, outdir)
}
