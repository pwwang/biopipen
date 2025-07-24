library(dplyr)

segfile = {{in.segfile | r}}
outfile = {{out.outfile | r}}
chrom_col = {{envs.chrom_col | r}}
excl_chroms = {{envs.excl_chroms | r}}
seg_col = {{envs.seg_col | r}}
segmean_transform = {{envs.segmean_transform | r}}

if (is.character(segmean_transform)) {
    segmean_transform = eval(parse(text=segmean_transform))
} # otherwise NULL


if (endsWith(segfile, ".vcf") || endsWith(segfile, ".vcf.gz")) {
  library(VariantAnnotation)
  segments = readVcf(segfile)
  seg = data.frame(
      chrom = as.character(seqnames(segments)),
      log2 = segments@info[[seg_col]]
  )
} else if (endsWith(segfile, ".bed")) {
  segments = read.table(segfile, header=F, row.names=NULL, sep="\t", stringsAsFactors=F)
  seg = data.frame(
      chrom = segments[, 1],
      log2 = segments[, 5]
  )
} else {
  segments = read.table(segfile, header=T, row.names=NULL, sep="\t", stringsAsFactors=F)
  seg = data.frame(
      chrom = segments[, chrom_col],
      log2 = segments[, seg_col]
  )
}
rm(segments)

if (!is.null(excl_chroms) && length(excl_chroms) > 0) {
    excl_chroms = sapply(excl_chroms, function(x) ifelse(startsWith(x, "chr"), x, paste0("chr", x)))
    seg = seg %>%
        mutate(.chrom = if_else(startsWith(chrom, "chr"), chrom, paste0("chr", chrom))) %>%
        filter(!.chrom %in% excl_chroms) %>%
        select(-.chrom)
}

if (!is.null(segmean_transform)) {
    seg$log2 = segmean_transform(seg$log2)
}

tmad = abs(mad(x = as.numeric(seg$log2), center=0, na.rm = TRUE))

cat(tmad, file=outfile)
