library(rlang)
library(rtracklayer)

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}
up <- {{envs.up | r}}
down <- {{envs.down | r}}
withbody <- {{envs.withbody | r}}
notfound <- {{envs.notfound | r}}
refgene <- {{envs.refgene | r}}
header <- {{envs.header | r}}
genecol <- {{envs.genecol | r}}
match_id <- {{envs.match_id | r}}
sort_ <- {{envs.sort | r}}
chrsize <- {{envs.chrsize | r}}

down <- down %||% up

data <- read.table(infile, header=header, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
genes <- data[[genecol]]
rm(data)

refgenes <- readGFF(refgene)
refcol <- ifelse(match_id, "gene_id", "gene_name")

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

# if withbody is FALSE:
#  make the promoters as start - up and start + down if strand is +
#  and end - down and end + up if strand is -
# if withbody is TRUE:
#  make the promoters as start - up and end + down if strand is +
#  and start - down and end + up if strand is -
promoters <- GRanges(
    seqnames = refgenes$seqid,
    ranges = IRanges(
        start = ifelse(refgenes$strand == "+", refgenes$start - up, refgenes$end - down),
        end = ifelse(refgenes$strand == "+", refgenes$start + down, refgenes$end + up)
    ),
    strand = refgenes$strand,
    name = refgenes[[refcol]]
)

if (withbody) {
    promoters$ranges <- IRanges(
        start = ifelse(refgenes$strand == "+", refgenes$start - up, refgenes$start - down),
        end = ifelse(refgenes$strand == "+", refgenes$end + down, refgenes$end + up)
    )
}

if (sort_) {
    chrom_sizes <- read.table(chrsize, header=FALSE, stringsAsFactors=FALSE, sep="\t")
    seqs <- chrom_sizes[chrom_sizes$V1 %in% seqlevels(promoters), 1, drop=TRUE]
    seqlevels(promoters) <- seqs
    promoters <- sort(promoters, ignore.strand = TRUE)
}

export.bed(promoters, outfile)
