library(AneuploidyScore)
library(dplyr)
library(tidyr)
library(tibble)
library(plotthis)
library(biopipen.utils)

segfile = {{in.segfile | r}}
outdir = {{out.outdir | r}}
chrom_col = {{envs.chrom_col | r}}
start_col = {{envs.start_col | r}}
end_col = {{envs.end_col | r}}
seg_col = {{envs.seg_col | r}}
cn_col = {{envs.cn_col | r}}
genome = {{envs.genome | r}}
threshold = {{envs.threshold | r}}
wgd_gf = {{envs.wgd_gf | r}}
excl_chroms = {{envs.excl_chroms | r}}

if (genome == "hg19") {
    data(ucsc.hg19.cytoband)
    cytoarm = cytobandToArm(ucsc.hg19.cytoband)
} else if (genome == "hg38") {
    data(ucsc.hg38.cytoband)
    cytoarm = cytobandToArm(ucsc.hg38.cytoband)
} else {
    stop(paste("Genome ", genome, " not supported"))
}

### Rewrite getCAA, as it raises error when some regions cannot be assigned
### to an arm
getCAA <- function(segf, cytoarm, tcn_col,
                   filter_centromere=FALSE, classifyCN=FALSE,
                   ploidy=0, threshold=0.2, ...){
  library(GenomicRanges)
  library(GenomeInfoDb)

  # tcn_col = 'Modal_Total_CN'
  stopifnot(AneuploidyScore:::.validateSeg(segf))

  ## Set up GenomicRange Objects of cytoband and seg files
  cyto_gr <- lapply(cytoarm, GenomicRanges::makeGRangesFromDataFrame,
                    keep.extra.columns=TRUE)
  seg_gr <- GenomicRanges::makeGRangesFromDataFrame(segf, keep.extra.columns = TRUE)
  seqlevelsStyle(seg_gr) <- 'UCSC'
  seg_chr <- split(seg_gr, f=seqnames(seg_gr))

  ## Goes through chr by chr,
  seg_cyto_chr <- lapply(names(seg_chr), function(chr_id){
    segc <- seg_chr[[chr_id]]
    cytoc <- sort(cyto_gr[[chr_id]])

    if(filter_centromere){
      ## Remove any segment that sligthly overlaps the centromere
      cen_ov <- findOverlaps(cytoc['cen',], segc)
      if(length(cen_ov) > 0) segc <- segc[-subjectHits(cen_ov),]
    }

    ## Create a GRanges object with all unique intervals between segc and cytoc
    starts <- tryCatch({
      sort(c(GenomicRanges::start(segc), GenomicRanges::start(cytoc)))
    }, error=function(e) {
      warning("Error to detect start on chromosome: ", chr_id, immediate. = TRUE)
      NULL
    })
    if (is.null(starts)) {
      return(NULL)
    }
    ends <- sort(c(GenomicRanges::end(segc), GenomicRanges::end(cytoc)))
    combc <- GRanges(seqnames=chr_id,
                     IRanges(start=unique(sort(c(starts, ends[-length(ends)]+1))),
                             end=unique(sort(c(ends, starts[-1]-1)))))

    # Map chr-arm to intervals
    cyto_comb_ov <- findOverlaps(cytoc, combc)
    ########### patched this line:
    combc <- combc[subjectHits(cyto_comb_ov),]
    mcols(combc)$arm <- names(cytoc[queryHits(cyto_comb_ov),])

    # Map seg values to intervals
    segc_comb_ov <- findOverlaps(segc, combc)
    mcols(combc)$CN <- NA
    mcols(combc[subjectHits(segc_comb_ov),])$CN <- mcols(segc[queryHits(segc_comb_ov),])[,tcn_col]

    ## Get chromosomal arm or whole-chromosome fractions for each CN interval
    .assembleFrac <- function(combc, assemble_method='arm', ...){
      if(assemble_method == 'arm'){
        X_fract <- lapply(split(combc, f=combc$arm), AneuploidyScore:::.getChrarmFractions, ...)
        ord <- order(factor(names(X_fract), levels=c("p", "cen", "q")))
        X_fract <- as.data.frame(do_call(rbind, X_fract[ord]))
        colnames(X_fract)[1:2] <- c("CAA_frac_NA", "CAA_frac_nonNA")
      } else {
        X_fract <- AneuploidyScore:::.getChrarmFractions(combc, ...)
        colnames(X_fract)[1:2] <- c("Chr_frac_NA", "Chr_frac_nonNA")
      }
      mcols(combc) <- cbind(mcols(combc), X_fract)
      return(combc)
    }

    combc <- .assembleFrac(combc, assemble_method='chr') # Chromosome fractions
    combc <- .assembleFrac(combc, assemble_method='arm', classifyCN=classifyCN,
                           ploidy=ploidy, threshold=threshold) # Chromosome arm fractions

    ## Handle intervals that have no CN value (NA; i.e. telomeric ends)
    combc$UID <- paste(combc$arm, round(combc$CN,2), sep="_")  ## Sets a unique arm_CN value
    if(any(is.na(combc$CN))) {
      combc_na <- combc[which(is.na(combc$CN)),]
      combc <- combc[-which(is.na(combc$CN)),] ## Removes intervals with no CN values (NA)
    }

    ## Reduce intervals where there is no CN change
    # e.g. c(4,4,4,2,4,4) => list(c(4,4,4), c(2), c(4,4))
    uid_int <- as.integer(factor(combc$UID))
    combc_arms <- split(combc, cumsum(c(TRUE, diff(uid_int) != 0)))
    combc_arms <- as(lapply(combc_arms, function(ca){
      ca_red <- reduce(ca)
      mcols(ca_red) <- mcols(ca)[1,]
      return(ca_red)
    }), "GRangesList")

    combc_arms <- unlist(combc_arms)
    mcols(combc_arms) <- mcols(combc_arms)[,-grep("^UID$", colnames(mcols(combc_arms)))]
    if(classifyCN){
      combc_arms$segCNclass <- AneuploidyScore:::.classifyCN(cn = combc_arms$CN, ploidy = ploidy,
                                           threshold=threshold)
    }


    return(combc_arms)
  })
  names(seg_cyto_chr) <- names(seg_chr)
  seg_cyto_chr <- seg_cyto_chr[!sapply(seg_cyto_chr, is.null)]
  return(as(seg_cyto_chr, "GRangesList"))
}

if (endsWith(segfile, ".vcf") || endsWith(segfile, ".vcf.gz")) {
  library(VariantAnnotation)
  vcf = readVcf(segfile)
  seg = data.frame(
      seqnames = as.character(seqnames(vcf)),
      start = start(vcf),
      end = vcf@info[[end_col]],
      seg.mean = vcf@info[[seg_col]]
  )
} else if (endsWith(segfile, ".bed")) {
  segments = read.table(segfile, header=F, row.names=NULL, sep="\t", stringsAsFactors=F)
  seg = data.frame(
      seqnames = segments[, 1],
      start = segments[, 2],
      end = segments[, 3],
      seg.mean = segments[, 5]
  )
} else {
  segments = read.table(segfile, header=T, row.names=NULL, sep="\t", stringsAsFactors=F)
  seg = data.frame(
      seqnames = segments[, chrom_col],
      start = segments[, start_col],
      end = segments[, end_col],
      seg.mean = segments[, seg_col]
  )
}

{% if envs.segmean_transform %}
segmean_transform = {{envs.segmean_transform}}
seg$seg.mean = segmean_transform(seg$seg.mean)
{% endif %}

if (!is.null(cn_col)) {
    seg$TCN = segments[, cn_col]
}

{% if envs.cn_transform %}
cn_transform = {{envs.cn_transform | r}}
if (is.character(cn_transform)) {
    cn_transform = eval(parse(text = cn_transform))
    seg$TCN = cn_transform(seg$seg.mean)
} else if (is.numeric(cn_transform) && length(cn_transform) > 1) {
    # Use cutoffs to transform
    # See also https://cnvkit.readthedocs.io/en/stable/pipeline.html#calling-methods
    # and https://github.com/etal/cnvkit/blob/9dd1e7c83705d1e1de6e6e4ab9fdc6973bf4002f/cnvlib/call.py#L98-L146
    # cn_transform =  c(-1.1, -0.25, 0.2, 0.7)
    # We assume neutral copy is 2
    # Still have to consider sex chromosomes
    neutral = 2
    TCN = cut(seg$seg.mean, breaks=c(-999, cn_transform), labels=seq(0, length(cn_transform) - 1))
    # Drop the factor structure
    TCN = as.integer(as.character(TCN))
    # Beyond the maximum cutoff, round it up
    TCN[is.na(TCN)] = ceiling(neutral * 2 ^ seg$seg.mean[is.na(TCN)])
    # If seg.mean is not available
    TCN[is.na(seg$seg.mean)] = neutral
    seg$TCN = TCN
}
{% endif %}

seg <- seg[
  !is.na(seg$seg.mean) & !is.na(seg$TCN) & !is.infinite(seg$seg.mean) & !is.infinite(seg$TCN),,
  drop=FALSE]

write.table(seg, file.path(outdir, "seg.txt"), sep="\t", quote=F, row.names=F, col.names=T)

wgd_ploidy = checkIfWGD(
    seg,
    tcn_col = "TCN",
    threshold = threshold,
    wgd_gf = wgd_gf,
    ploidy_method = "wmean"
)

seg_caa = getCAA(
    seg,
    cytoarm,
    tcn_col = "TCN",
    classifyCN = TRUE,
    ploidy = wgd_ploidy['ploidy'],
    threshold = threshold
)

norm_chrom = function(chr) {
    if (grepl("^chr", chr)) {
        chr = gsub("^chr", "", chr)
    }
    return(chr)
}

caa = reduceArms(seg_caa, caa_method = c("arm", "seg"), arm_ids=c("p", "q"))
if (length(excl_chroms) > 0) {
  excl_chroms = paste0(excl_chroms, "_")
  excluded = sapply(rownames(caa), function(chr) {
    any(sapply(excl_chroms, function(excl_chr) {
      startsWith(norm_chrom(chr), norm_chrom(excl_chr))
    }))
  })
  AS = colSums(abs(caa[!excluded, ]), na.rm = TRUE)
} else {
  AS = colSums(abs(caa), na.rm = TRUE)
}
caa = caa %>% as.data.frame() %>% tibble::rownames_to_column("Arms")

write.table(caa, file.path(outdir, "CAA.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(as.data.frame(AS), file.path(outdir, "AS.txt"), sep="\t", quote=F, row.names=T, col.names=F)

# plots
plotdata = caa |>
    extract(Arms, c("arm_id", "arm_name"), "chr(\\d+)_(p|q)", remove=FALSE) |>
    mutate(arm_id = as.numeric(arm_id)) |>
    arrange(arm_id, arm_name) |>
    mutate(
        Type=case_when(seg > 0 ~ "Duplication", seg < 0 ~ "Deletion", TRUE ~ "Neutral"),
        Arms=factor(Arms, levels=Arms)
    ) |>
    select(Arms, Type, arm, seg) |>
    pivot_longer(c("arm", "seg"), names_to="SignalType", values_to="Signal")

sig_min = min(-1, plotdata$Signal, na.rm=TRUE)
sig_max = max(1, plotdata$Signal, na.rm=TRUE)

png(file.path(outdir, "AneuploidyScore.png"), width=1000, height=600, res=100)
p <- BarPlot(
    plotdata,
    x = "Arms",
    y = "Signal",
    fill = "Type",
    facet_by = "SignalType",
    facet_nrow = 2,
    y_min = sig_min,
    y_max = sig_max,
    x_text_angle = 90,
    aspect.ratio = 0.2
)
print(p)
dev.off()
