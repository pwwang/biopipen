library(AneuploidyScore)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggprism)

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

if (genome == "hg19") {
    data(ucsc.hg19.cytoband)
    cytoarm = cytobandToArm(ucsc.hg19.cytoband)
} else if (genome == "hg38") {
    data(ucsc.hg38.cytoband)
    cytoarm = cytobandToArm(ucsc.hg38.cytoband)
} else {
    stop(paste("Genome ", genome, " not supported"))
}

{% if envs.segmean_transform %}
segmean_transform = {{envs.segmean_transform}}
{% else %}
segmean_transform = NULL
{% endif %}

{% if envs.cn_transform %}
cn_transform = {{envs.cn_transform}}
{% else %}
cn_transform = NULL
{% endif %}

segments = read.table(segfile, header=T, row.names=NULL, sep="\t", stringsAsFactors=F)
seg = data.frame(
    seqnames = segments[, chrom_col],
    start = segments[, start_col],
    end = segments[, end_col],
    seg.mean = segments[, seg_col]
)

if (!is.null(cn_transform) && !is.null(cn_col)) {
    seg$TCN = cn_transform(seg$seg.mean)
}

if (!is.null(cn_col)) {
    seg$TCN = segments[, cn_col]
}

if (!is.null(segmean_transform)) {
    seg$seg.mean = segmean_transform(seg$seg.mean)
}

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

caa = reduceArms(seg_caa, caa_method = c("arm", "seg"), arm_ids=c("p", "q"))
AS = colSums(abs(caa), na.rm = TRUE)
caa = caa %>% as.data.frame() %>% tibble::rownames_to_column("Arms")

write.table(caa, file.path(outdir, "CAA.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(as.data.frame(AS), file.path(outdir, "AS.txt"), sep="\t", quote=F, row.names=T, col.names=F)

# plots
all_arms = c(
    paste0("chr", c(1:22, "X", "Y"), "_p"),
    paste0("chr", c(1:22, "X", "Y"), "_q")
)
plotdata = caa
missing_arms = setdiff(all_arms, plotdata$Arms)
for (arm in missing_arms) {
    plotdata = bind_rows(plotdata, list(Arms=arm, arm=NA, seg=NA))
}

plotdata = plotdata |>
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
ggplot(plotdata) +
    geom_bar(aes(x=Arms, y=Signal, fill=Type), stat="identity") +
    geom_hline(yintercept=0, color="black", size=0.1) +
    ylim(c(sig_min, sig_max)) +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~SignalType, scales="free_y", nrow=2)
dev.off()
