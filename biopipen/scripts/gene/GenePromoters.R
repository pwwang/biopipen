library(rlang)
library(rtracklayer)

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}
up <- {{envs.up | r}}
down <- {{envs.down | r}}
notfound <- {{envs.notfound | r}}
refgene <- {{envs.refgene | r}}
header <- {{envs.header | r}}
genecol <- {{envs.genecol | r}}
match_id <- {{envs.match_id | r}}
sort_ <- {{envs.sort | r}}
chrsize <- {{envs.chrsize | r}}

down <- down %||% up

refgenes <- readGFF(refgene)
refcol <- ifelse(match_id, "gene_id", "gene_name")

if (infile == "/dev/null") {
    genes <- unique(refgenes[[refcol]])
} else {
    data <- read.table(infile, header=header, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
    genes <- data[[genecol]]
    rm(data)
}

notfound_genes <- setdiff(genes, refgenes[[refcol]])
if (notfound == "error" && length(notfound_genes) > 0) {
    stop(paste(
        "The following genes were not found in the reference annotation:",
        paste(notfound_genes, collapse=", ")
    ))
} else if (notfound == 'skip') {
    genes <- genes[!genes %in% notfound_genes]
}

# Select the genes that are in the reference annotation and keep the order
# of the records in genes
refgenes <- refgenes[match(genes, refgenes[[refcol]]), , drop = FALSE]
refgenes <- unique(makeGRangesFromDataFrame(refgenes, keep.extra.columns=TRUE))

proms <- promoters(refgenes, up=up, down=down)
# Scores must be non-NA numeric values
elementMetadata(proms)$name <- elementMetadata(proms)[[refcol]]
score(proms) <- 0
start(proms) <- pmax(1, start(proms))

if (sort_) {
    chrom_sizes <- read.table(chrsize, header=FALSE, stringsAsFactors=FALSE, sep="\t")
    common_chroms <- intersect(chrom_sizes$V1, seqlevels(proms))
    if (length(common_chroms) == 0) {
        stop("No common chromosomes found between the promoters and the chromosome sizes. Do you use the correct chromosome sizes file?")
    }
    proms <- keepSeqlevels(proms, common_chroms, pruning.mode="coarse")
    seqlevels(proms) <- common_chroms
    proms <- sort(proms, ignore.strand = TRUE)
}

export.bed(proms, outfile)
