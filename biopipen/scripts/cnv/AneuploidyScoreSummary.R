library(dplyr)
library(tidyr)
library(tibble)
library(plotthis)
library(biopipen.utils)

asdirs <- {{in.asdirs | r}}
metafile <- {{in.metafile | r}}
outdir <- {{out.outdir | r}}
group_cols <- {{envs.group_cols | r}}
sample_name_fun <- {{envs.sample_name | r}}
heatmap_cases <- {{envs.heatmap_cases | r}}

if (!is.null(sample_name_fun)) {
    sample_name_fun <- eval(parse(text = sample_name_fun))
}

get_sample_from_asdir <- function(asdir) {
    x <- basename(asdir)
    if (endsWith(x, ".aneuploidy_score")) {
        x <- substr(x, 1, nchar(x) - 17)
    }
    if (endsWith(x, ".call")) {
        x <- substr(x, 1, nchar(x) - 5)
    }
    if (!is.null(sample_name_fun)) {
        x <- sample_name_fun(x)
    }
    x
}

asdir_to_sample <- lapply(asdirs, get_sample_from_asdir)
names(asdir_to_sample) <- asdirs
table_sams <- table(unlist(asdir_to_sample))
if (any(table_sams > 1)) {
    log_warn("Duplicate sample names found in asdirs: ")
    dup_sams <- names(table_sams[table_sams > 1])
    for (dup_sam in dup_sams) {
        i <- 1
        for (asdir in asdirs) {
            if (asdir_to_sample[[asdir]] == dup_sam) {
                dedup_sam <- paste0(dup_sam, "_", i)
                log_warn(paste0("- Changing ", dup_sam, "(", asdir, ") to ", dedup_sam))
                asdir_to_sample[[asdir]] <- paste0(dup_sam, "_", i)
                i <- i + 1
            }
        }
    }
}
sams <- unlist(asdir_to_sample)

meta_cols <- c()
if (!is.null(group_cols)) {
    for (group_col in group_cols) {
        if (grepl(",", group_col, fixed = TRUE)) {
            subcols <- strsplit(group_col, ",")[[1]]
            if (length(subcols) > 2) {
                stop("Only support 2 columns combined for group_cols")
            }
            meta_cols <- union(meta_cols, subcols)
        } else {
            meta_cols <- union(meta_cols, group_col)
        }
    }
}

if (!is.null(metafile)) {
    metadf <- read.table(metafile, header=T, row.names=NULL, sep="\t", stringsAsFactors=F)
    if (!is.null(metadf$Sample)) {
        metadf$Sample <- as.character(metadf$Sample)
    } else {
        colnames(metadf)[1] <- "Sample"
    }
    metadf <- metadf[metadf$Sample %in% sams, c("Sample", meta_cols), drop=FALSE]
    rownames(metadf) <- metadf$Sample
    if (nrow(metadf) != length(sams)) {
        stop(paste("Not all samples in metafile:", paste(setdiff(sams, metadf$Sample), collapse=", ")))
    }
} else {
    metadf <- NULL
    if (!is.null(group_cols) && length(group_cols) > 0) {
        stop("`envs.group_cols` given but no metafile provided")
    }
}

read_caa <- function(asdir) {
    # Sample Arms arm seg
    sample <- asdir_to_sample[[asdir]]
    caa <- read.table(
        file.path(asdir, "CAA.txt"),
        header=T,
        row.names=NULL,
        sep="\t",
        stringsAsFactors=F,
    )
    caa$Sample <- sample
    caa
}

read_as <- function(asdir) {
    # Sample SignalType Signal
    sample <- asdir_to_sample[[asdir]]
    as <- read.table(
        file.path(asdir, "AS.txt"),
        header=F,
        row.names=NULL,
        sep="\t",
        stringsAsFactors=F,
    )
    colnames(as) <- c("SignalType", "Signal")
    as$Sample <- sample
    as
}

# Sample Arms arm seg
caa <- do_call(rbind, lapply(asdirs, read_caa))
# Sample SignalType Signal
as <- do_call(rbind, lapply(asdirs, read_as))

# Sample chr1_p chr1_q chr2_p chr2_q ...
caa_arm <- caa %>%
    select(-"seg") %>%
    pivot_wider(names_from="Arms", values_from="arm")

# Sample chr1_p chr1_q chr2_p chr2_q ...
caa_seg <- caa %>%
    select(-"arm") %>%
    pivot_wider(names_from="Arms", values_from="seg")

# Sample SignalType Signal
as_arm <- as %>% filter(SignalType == "arm") %>% select(-"SignalType")
as_seg <- as %>% filter(SignalType == "seg") %>% select(-"SignalType")

if (!is.null(metadf)) {
    caa_arm <- caa_arm %>% left_join(metadf, by="Sample")
    caa_seg <- caa_seg %>% left_join(metadf, by="Sample")
    as_arm <- as_arm %>% left_join(metadf, by="Sample")
    as_seg <- as_seg %>% left_join(metadf, by="Sample")
}

write.table(caa_arm, file.path(outdir, "CAA_arm.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(caa_seg, file.path(outdir, "CAA_seg.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(as_arm, file.path(outdir, "AS_arm.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(as_seg, file.path(outdir, "AS_seg.txt"), sep="\t", quote=F, row.names=F, col.names=T)

# Plot AS without grouping
p_as_arm <- BarPlot(
    as_arm,
    x="Sample",
    y="Signal",
    title="Aneuploidy Score (Arm)",
    xlab="Sample",
    ylab="Aneuploidy Score",
    x_text_angle = 90
)

png(file.path(outdir, "AS_arm.png"), width=400 + nrow(caa_arm) * 12, height=600, res=100)
print(p_as_arm)
dev.off()

p_as_seg <- BarPlot(
    as_seg,
    x="Sample",
    y="Signal",
    title="Aneuploidy Score (Segment)",
    xlab="Sample",
    ylab="Aneuploidy Score",
    x_text_angle = 90
)

png(file.path(outdir, "AS_seg.png"), width=400 + nrow(caa_seg) * 12, height=600, res=100)
print(p_as_seg)
dev.off()

# Plot AS for each group_col
if (!is.null(group_cols)) {
    for (group_col in group_cols) {
        if (!grepl(",", group_col, fixed = TRUE)) {

            p_as_arm_bar_group <- BarPlot(
                as_arm,
                x="Sample",
                y="Signal",
                fill=group_col,
                title=paste0("Aneuploidy Score (Arm) - ", group_col),
                xlab="Sample",
                ylab="Aneuploidy Score",
                x_text_angle = 90
            )

            png(file.path(outdir, paste0("AS_arm_bar_", group_col, ".png")), width=400 + nrow(caa_arm) * 12, height=600, res=100)
            print(p_as_arm_bar_group)
            dev.off()

            p_as_seg_bar_group <- BarPlot(
                as_seg,
                x="Sample",
                y="Signal",
                fill=group_col,
                title=paste0("Aneuploidy Score (Segment) - ", group_col),
                xlab="Sample",
                ylab="Aneuploidy Score",
                x_text_angle = 90
            )

            png(file.path(outdir, paste0("AS_seg_bar_", group_col, ".png")), width=400 + nrow(caa_seg) * 12, height=600, res=100)
            print(p_as_seg_bar_group)
            dev.off()

            # Voilin + boxplot

            p_as_arm_violin_group <- ViolinPlot(
                as_arm,
                x=group_col,
                y="Signal",
                title=paste0("Aneuploidy Score (Arm) - ", group_col),
                xlab=group_col,
                ylab="Aneuploidy Score",
                x_text_angle = 90,
                comparisons = TRUE,
                sig_label = "p.format",
                add_point = TRUE,
                add_box = TRUE
            )

            png(file.path(outdir, paste0("AS_arm_violin_", group_col, ".png")), width=1000, height=600, res=100)
            print(p_as_arm_violin_group)
            dev.off()

            p_as_seg_violin_group <- ViolinPlot(
                as_seg,
                x=group_col,
                y="Signal",
                title=paste0("Aneuploidy Score (Segment) - ", group_col),
                xlab=group_col,
                ylab="Aneuploidy Score",
                x_text_angle = 90,
                comparisons = TRUE,
                sig_label = "p.format",
                add_point = TRUE,
                add_box = TRUE
            )

            png(file.path(outdir, paste0("AS_seg_violin_", group_col, ".png")), width=1000, height=600, res=100)
            print(p_as_seg_violin_group)
            dev.off()

        } else {
            # Multiple layers with group_col
            group_cols <- strsplit(group_col, ",")[[1]]
            group_col1 <- group_cols[1]
            group_col2 <- group_cols[2]

            # For each group_col1, plot a barplot with group_col2 as fill, and
            # concatenate them together using patch work, with ncol=2
            # calcuate the height and width of the plot based on the number of
            # groups
            as_arm <- as_arm %>% arrange(!!sym(group_col1), !!sym(group_col2)) %>% mutate(Sample=factor(Sample, levels=Sample))
            p <- BarPlot(
                as_arm,
                x="Sample",
                y="Signal",
                split_by=group_col1,
                fill=group_col2,
                xlab="Sample",
                ylab="Aneuploidy Score",
                x_text_angle = 90,
                ncol = 2
            )

            png(
                file.path(outdir, paste0("AS_arm_bar_", group_col, ".png")),
                width=1000,
                height=length(unique(as_arm[[group_col1]])) * 200,
                res=100
            )
            print(p)
            dev.off()

            as_seg <- as_seg %>% arrange(!!sym(group_col1), !!sym(group_col2)) %>% mutate(Sample=factor(Sample, levels=Sample))
            p <- BarPlot(
                as_seg,
                x="Sample",
                y="Signal",
                split_by=group_col1,
                fill=group_col2,
                xlab="Sample",
                ylab="Aneuploidy Score",
                x_text_angle = 90,
                ncol = 2
            )
            png(
                file.path(outdir, paste0("AS_seg_bar_", group_col, ".png")),
                width=1000,
                height=length(unique(as_seg[[group_col1]])) * 200,
                res=100
            )
            print(p)
            dev.off()

            # Do the same for Voilin + boxplot
            p <- ViolinPlot(
                as_arm,
                x=group_col2,
                y="Signal",
                split_by = group_col1,
                xlab=group_col2,
                ylab="Aneuploidy Score",
                x_text_angle = 90,
                comparisons = TRUE,
                sig_label = "p.format",
                add_point = TRUE,
                add_box = TRUE,
                ncol = 2
            )

            png(
                file.path(outdir, paste0("AS_arm_violin_", group_col, ".png")),
                width=1000,
                height=length(unique(as_arm[[group_col1]])) * 200,
                res=100
            )
            print(p)
            dev.off()

            p <- ViolinPlot(
                as_seg,
                x=group_col2,
                y="Signal",
                split_by = group_col1,
                xlab=group_col2,
                ylab="Aneuploidy Score",
                x_text_angle = 90,
                comparisons = TRUE,
                sig_label = "p.format",
                add_point = TRUE,
                add_box = TRUE,
                ncol = 2
            )

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
    arms <- heatmap_cases[[heatmap_name]]
    if (all(arms != "ALL")) {
        caa_df <- caa_arm %>% select(Sample, !!meta_cols, !!arms)
    } else {
        caa_df <- caa_arm
    }
    caa_df <- caa_df %>% column_to_rownames("Sample")
    if (!is.null(metadf)) {
        caa_df <- caa_df %>% select(-!!meta_cols)
    }
    caa_df <- caa_df %>%
        rownames_to_column("Sample") %>%
        pivot_longer(cols=-"Sample", names_to="Arms", values_to="Signal") %>%
        pivot_wider(names_from="Sample", values_from="Signal")

    height <- 300 + 20 * ncol(caa_df)  # all arms: 300 + 30 * 46 = 1680
    width <- 300 + 25 * nrow(caa_df)  # 10 samples: 300 + 30 * 10 = 600
    # print(caa_df)
    hmp <- Heatmap(
        caa_df,
        rows_data = metadf,
        name = "CAA",
        rows_by = setdiff(colnames(caa_df), "Arms"),
        columns_by = "Arms",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_annotation = colnames(metadf),
        lower_cutoff = -1,
        upper_cutoff = 1
    )

    png(
        file.path(outdir, paste0("Heatmap_", heatmap_name, ".png")),
        width=width,
        height=height,
        res=100
    )
    plot(hmp)
    dev.off()
}
