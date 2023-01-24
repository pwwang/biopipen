library(ggplot2)
library(ggprism)
library(dplyr)
library(tidyr)
library(tibble)
source("{{biopipen_dir}}/utils/plot.R")

asdirs = {{in.asdirs | r}}
metafile  = {{in.metafile | r}}
outdir = {{out.outdir | r}}
group_col = {{envs.group_col | r}}
heatmap_cases = {{envs.heatmap_cases | r}}

meta = NULL
if (file.exists(metafile)) {
    metadf = read.table(metafile, header=T, row.names=NULL, sep="\t", stringsAsFactors=F)
    meta = as.list(metadf[[group_col]])
    names(meta) = metadf[, 1, drop=TRUE]
}

stem0 = function(path) {
    strsplit(basename(path), ".", fixed=TRUE)[[1]][1]
}

read_caa = function(asdir) {
    sample = stem0(asdir)
    caa = read.table(
        file.path(asdir, "CAA.txt"),
        header=T,
        row.names=NULL,
        sep="\t",
        stringsAsFactors=F,
    )
    caa$Sample = sample
    if (!is.null(meta)) {
        caa$Group = meta[[sample]]
    }
    caa
}

read_as = function(asdir) {
    sample = stem0(asdir)
    as = read.table(
        file.path(asdir, "AS.txt"),
        header=F,
        row.names=NULL,
        sep="\t",
        stringsAsFactors=F,
    )
    colnames(as) = c("SignalType", "Signal")
    as$Sample = sample
    if (!is.null(meta)) {
        as$Group = meta[[sample]]
    }
    as
}

caa = do.call(rbind, lapply(asdirs, read_caa))
as = do.call(rbind, lapply(asdirs, read_as))

caa_arm = caa %>%
    select(-"seg") %>%
    pivot_wider(names_from="Arms", values_from="arm")

caa_seg = caa %>%
    select(-"arm") %>%
    pivot_wider(names_from="Arms", values_from="seg")

if (!is.null(meta)) {
    caa_arm = caa_arm %>% arrange(Group)
    caa_seg = caa_seg %>% arrange(Group)
    plotdata = as %>%
        arrange(Group) %>%
        mutate(Sample = factor(Sample, levels=unique(Sample)))
} else {
    plotdata = as
}

write.table(caa_arm, file.path(outdir, "CAA_arm.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(caa_seg, file.path(outdir, "CAA_seg.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(plotdata %>% filter(SignalType == "arm"), file.path(outdir, "AS_arm.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(plotdata %>% filter(SignalType == "seg"), file.path(outdir, "AS_seg.txt"), sep="\t", quote=F, row.names=F, col.names=T)

png(file.path(outdir, "CAAs.png"), width=1000, height=600, res=100)
if (!is.null(meta)) {
    mapping = aes_string(x="Sample", y="Signal", fill="Group")
} else {
    mapping = aes_string(x="Sample", y="Signal")
}
ggplot(plotdata) +
    geom_bar(mapping, stat="identity") +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~SignalType, scales="free_y", nrow=2)
dev.off()


if (!is.null(meta)) {
    png(file.path(outdir, "CAAs_group.png"), width=800, height=800, res=100)
    p = ggplot(plotdata) +
        geom_violin(aes(x=Group, y=Signal, fill=Group)) +
        geom_boxplot(aes(x=Group, y=Signal), width=.1) +
        theme_prism() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        facet_wrap(~SignalType, scales="free_y", nrow=2)
    print(p)
    dev.off()
}

# Heatmaps
for (heatmap_name in names(heatmap_cases)) {
    arms = heatmap_cases[[heatmap_name]]
    if (all(arms != "ALL")) {
        caa_df = caa_arm %>% select(Sample, Group, !!arms)
    } else {
        caa_df = caa_arm
    }
    caa_df = caa_df %>% column_to_rownames("Sample")
    if (!is.null(meta)) {
        caa_df = caa_df %>% select(-"Group")
    }

    width = 300 + 20 * ncol(caa_df)  # all arms: 300 + 30 * 46 = 1680
    height = 300 + 30 * nrow(caa_df)  # 10 samples: 300 + 30 * 10 = 600
    args = list(
        name = "CAA",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_names_side = "left",
        rect_gp = grid::gpar(col = "#FFFFFF", lwd = 1)
    )
    if (!is.null(meta)) {
        args$right_annotation = ComplexHeatmap::rowAnnotation(Group = caa_arm$Group)
    }
    plotHeatmap(
        caa_df,
        args = args,
        devpars = list(width=width, height=height, res=100),
        outfile = file.path(outdir, paste0("Heatmap_", heatmap_name, ".png"))
    )
}
