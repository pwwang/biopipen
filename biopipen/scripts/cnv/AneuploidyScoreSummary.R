source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/plot.R")

library(ggplot2)
library(ggprism)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)

asdirs = {{in.asdirs | r}}
metafile  = {{in.metafile | r}}
outdir = {{out.outdir | r}}
group_cols = {{envs.group_cols | r}}
sample_name_fun = {{envs.sample_name | r}}
heatmap_cases = {{envs.heatmap_cases | r}}

if (!is.null(sample_name_fun)) {
    sample_name_fun = eval(parse(text=sample_name_fun))
}

get_sample_from_asdir = function(asdir) {
    x = basename(asdir)
    if (endsWith(x, ".aneuploidy_score")) {
        x = substr(x, 1, nchar(x) - 17)
    }
    if (endsWith(x, ".call")) {
        x = substr(x, 1, nchar(x) - 5)
    }
    if (!is.null(sample_name_fun)) {
        x = sample_name_fun(x)
    }
    x
}

sams = sapply(asdirs, get_sample_from_asdir)

meta_cols = c()
if (!is.null(group_cols)) {
    for (group_col in group_cols) {
        if (grepl(",", group_col, fixed = TRUE)) {
            subcols = strsplit(group_col, ",")[[1]]
            if (length(subcols) > 2) {
                stop("Only support 2 columns combined for group_cols")
            }
            meta_cols = union(meta_cols, subcols)
        } else {
            meta_cols = union(meta_cols, group_col)
        }
    }
}

if (!is.null(metafile)) {
    metadf = read.table(metafile, header=T, row.names=NULL, sep="\t", stringsAsFactors=F)
    if (!is.null(metadf$Sample)) {
        metadf$Sample = as.character(metadf$Sample)
    } else {
        colnames(metadf)[1] = "Sample"
    }
    metadf = metadf[metadf$Sample %in% sams, c("Sample", meta_cols), drop=FALSE]
    if (nrow(metadf) != length(sams)) {
        stop(paste("Not all samples in metafile:", paste(setdiff(sams, metadf$Sample), collapse=", ")))
    }
} else {
    metadf = NULL
    if (!is.null(group_cols) && length(group_cols) > 0) {
        stop("`envs.group_cols` given but no metafile provided")
    }
}



read_caa = function(asdir) {
    # Sample Arms arm seg
    sample = get_sample_from_asdir(asdir)
    caa = read.table(
        file.path(asdir, "CAA.txt"),
        header=T,
        row.names=NULL,
        sep="\t",
        stringsAsFactors=F,
    )
    caa$Sample = sample
    caa
}

read_as = function(asdir) {
    # Sample SignalType Signal
    sample = get_sample_from_asdir(asdir)
    as = read.table(
        file.path(asdir, "AS.txt"),
        header=F,
        row.names=NULL,
        sep="\t",
        stringsAsFactors=F,
    )
    colnames(as) = c("SignalType", "Signal")
    as$Sample = sample
    as
}

# Sample Arms arm seg
caa = do_call(rbind, lapply(asdirs, read_caa))
# Sample SignalType Signal
as = do_call(rbind, lapply(asdirs, read_as))

# Sample chr1_p chr1_q chr2_p chr2_q ...
caa_arm = caa %>%
    select(-"seg") %>%
    pivot_wider(names_from="Arms", values_from="arm")

# Sample chr1_p chr1_q chr2_p chr2_q ...
caa_seg = caa %>%
    select(-"arm") %>%
    pivot_wider(names_from="Arms", values_from="seg")

# Sample SignalType Signal
as_arm = as %>% filter(SignalType == "arm") %>% select(-"SignalType")
as_seg = as %>% filter(SignalType == "seg") %>% select(-"SignalType")

if (!is.null(metadf)) {
    caa_arm = caa_arm %>% left_join(metadf, by="Sample")
    caa_seg = caa_seg %>% left_join(metadf, by="Sample")
    as_arm = as_arm %>% left_join(metadf, by="Sample")
    as_seg = as_seg %>% left_join(metadf, by="Sample")
}


write.table(caa_arm, file.path(outdir, "CAA_arm.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(caa_seg, file.path(outdir, "CAA_seg.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(as_arm, file.path(outdir, "AS_arm.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(as_seg, file.path(outdir, "AS_seg.txt"), sep="\t", quote=F, row.names=F, col.names=T)

# Plot AS without grouping
p_as_arm = ggplot(as_arm) +
    geom_bar(aes(x=Sample, y=Signal), stat="identity") +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

png(file.path(outdir, "AS_arm.png"), width=400 + nrow(caa_arm) * 12, height=600, res=100)
print(p_as_arm)
dev.off()

p_as_seg = ggplot(as_seg) +
    geom_bar(aes(x=Sample, y=Signal), stat="identity") +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

png(file.path(outdir, "AS_seg.png"), width=400 + nrow(caa_seg) * 12, height=600, res=100)
print(p_as_seg)
dev.off()

# Plot AS for each group_col
if (!is.null(group_cols)) {
    for (group_col in group_cols) {
        if (!grepl(",", group_col, fixed = TRUE)) {
            # Single layer with group_col
            p_as_arm_bar_group = ggplot(
                    as_arm %>% arrange(!!sym(group_col)) %>% mutate(Sample=factor(Sample, levels=Sample))
                ) +
                geom_bar(aes(x=Sample, y=Signal, fill=!!sym(group_col)), stat="identity") +
                theme_prism() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1))

            png(file.path(outdir, paste0("AS_arm_bar_", group_col, ".png")), width=400 + nrow(caa_arm) * 12, height=600, res=100)
            print(p_as_arm_bar_group)
            dev.off()

            p_as_seg_bar_group = ggplot(
                    as_seg %>% arrange(!!sym(group_col)) %>% mutate(Sample=factor(Sample, levels=Sample))
                ) +
                geom_bar(aes(x=Sample, y=Signal, fill=!!sym(group_col)), stat="identity") +
                theme_prism() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1))

            png(file.path(outdir, paste0("AS_seg_bar_", group_col, ".png")), width=400 + nrow(caa_seg) * 12, height=600, res=100)
            print(p_as_seg_bar_group)
            dev.off()

            # Voilin + boxplot
            p_as_arm_violin_group = ggplot(
                    as_arm %>% arrange(!!sym(group_col)) %>% mutate(Sample=factor(Sample, levels=Sample))
                ) +
                geom_violin(aes(x=!!sym(group_col), y=Signal), fill="steelblue", trim=FALSE) +
                geom_boxplot(aes(x=!!sym(group_col), y=Signal), width=0.1, outlier.shape=NA) +
                theme_prism() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1))

            png(file.path(outdir, paste0("AS_arm_violin_", group_col, ".png")), width=1000, height=600, res=100)
            print(p_as_arm_violin_group)
            dev.off()

            p_as_seg_violin_group = ggplot(
                    as_seg %>% arrange(!!sym(group_col)) %>% mutate(Sample=factor(Sample, levels=Sample))
                ) +
                geom_violin(aes(x=!!sym(group_col), y=Signal), fill="steelblue", trim=FALSE) +
                geom_boxplot(aes(x=!!sym(group_col), y=Signal), width=0.1, outlier.shape=NA) +
                theme_prism() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1))

            png(file.path(outdir, paste0("AS_seg_violin_", group_col, ".png")), width=1000, height=600, res=100)
            print(p_as_seg_violin_group)
            dev.off()

        } else {
            # Multiple layers with group_col
            group_cols = strsplit(group_col, ",")[[1]]
            group_col1 = group_cols[1]
            group_col2 = group_cols[2]

            # For each group_col1, plot a barplot with group_col2 as fill, and
            # concatenate them together using patch work, with ncol=2
            # calcuate the height and width of the plot based on the number of
            # groups
            ps = as_arm %>%
                group_by(!!sym(group_col1)) %>%
                group_map(function(.x, .y) {
                    p = ggplot(.x %>% arrange(!!sym(group_col2)) %>% mutate(Sample=factor(Sample, levels=Sample))) +
                        geom_bar(aes(x=Sample, y=Signal, fill=!!sym(group_col2)), stat="identity") +
                        theme_prism() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                        ggtitle(.y[[group_col1]][1])
                    p
                })

            p = wrap_plots(ps, ncol=2)
            png(
                file.path(outdir, paste0("AS_arm_bar_", group_col, ".png")),
                width=1000,
                height=length(unique(as_arm[[group_col1]])) * 200,
                res=100
            )
            print(p)
            dev.off()

            ps = as_seg %>%
                group_by(!!sym(group_col1)) %>%
                group_map(function(.x, .y) {
                    p = ggplot(.x %>% arrange(!!sym(group_col2)) %>% mutate(Sample=factor(Sample, levels=Sample))) +
                        geom_bar(aes(x=Sample, y=Signal, fill=!!sym(group_col2)), stat="identity") +
                        theme_prism() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                        ggtitle(.y[[group_col1]][1])
                    p
                })

            p = wrap_plots(ps, ncol=2)
            png(
                file.path(outdir, paste0("AS_seg_bar_", group_col, ".png")),
                width=1000,
                height=length(unique(as_seg[[group_col1]])) * 200,
                res=100
            )
            print(p)
            dev.off()

            # Do the same for Voilin + boxplot
            ps = as_arm %>%
                group_by(!!sym(group_col1)) %>%
                group_map(function(.x, .y) {
                    p = ggplot(.x %>% arrange(!!sym(group_col2)) %>% mutate(Sample=factor(Sample, levels=Sample))) +
                        geom_violin(aes(x=!!sym(group_col2), y=Signal), fill="steelblue", trim=FALSE) +
                        geom_boxplot(aes(x=!!sym(group_col2), y=Signal), width=0.1, outlier.shape=NA) +
                        theme_prism() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                        ggtitle(.y[[group_col1]][1])
                    p
                })

            p = wrap_plots(ps, ncol=2)
            png(
                file.path(outdir, paste0("AS_arm_violin_", group_col, ".png")),
                width=1000,
                height=length(unique(as_arm[[group_col1]])) * 200,
                res=100
            )
            print(p)
            dev.off()

            ps = as_seg %>%
                group_by(!!sym(group_col1)) %>%
                group_map(function(.x, .y) {
                    p = ggplot(.x %>% arrange(!!sym(group_col2)) %>% mutate(Sample=factor(Sample, levels=Sample))) +
                        geom_violin(aes(x=!!sym(group_col2), y=Signal), fill="steelblue", trim=FALSE) +
                        geom_boxplot(aes(x=!!sym(group_col2), y=Signal), width=0.1, outlier.shape=NA) +
                        theme_prism() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                        ggtitle(.y[[group_col1]][1])
                    p
                })

            p = wrap_plots(ps, ncol=2)
            png(
                file.path(outdir, paste0("AS_seg_violin_", group_col, ".png")),
                width=1000,
                height=length(unique(as_seg[[group_col1]])) * 200,
                res=100
            )
            print(p)
            dev.off()
        }
    }
}

# Heatmaps
for (heatmap_name in names(heatmap_cases)) {
    arms = heatmap_cases[[heatmap_name]]
    if (all(arms != "ALL")) {
        caa_df = caa_arm %>% select(Sample, !!meta_cols, !!arms)
    } else {
        caa_df = caa_arm
    }
    caa_df = caa_df %>% column_to_rownames("Sample")
    if (!is.null(metadf)) {
        caa_df = caa_df %>% select(-!!meta_cols)
    }

    width = 300 + 20 * ncol(caa_df)  # all arms: 300 + 30 * 46 = 1680
    height = 300 + 25 * nrow(caa_df)  # 10 samples: 300 + 30 * 10 = 600
    args = list(
        name = "CAA",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_names_side = "left",
        rect_gp = grid::gpar(col = "#FFFFFF", lwd = 1)
    )
    if (!is.null(metadf)) {
        row_annos = list()
        for (meta_col in meta_cols) {
            row_annos[[meta_col]] = metadf[[meta_col]]
        }
        if (length(row_annos) > 0) {
            args$right_annotation = do_call(ComplexHeatmap::rowAnnotation, row_annos)
        }
    }
    plotHeatmap(
        caa_df,
        args = args,
        devpars = list(width=width, height=height, res=100),
        outfile = file.path(outdir, paste0("Heatmap_", heatmap_name, ".png"))
    )
}
