library(dplyr)
library(tidyr)
library(tibble)
library(plotthis)

tmadfiles = {{in.tmadfiles | r}}
metafile  = {{in.metafile | r}}
outdir = {{out.outdir | r}}
group_cols = {{envs.group_cols | r}}
sample_name_fun = {{envs.sample_name | r}}

if (!is.null(sample_name_fun)) {
    sample_name_fun = eval(parse(text=sample_name_fun))
}

read_tmad = function(tmadfile) {
    as.numeric(suppressWarnings(readLines(tmadfile)))
}

tmads = sapply(tmadfiles, read_tmad)
sams = sapply(tmadfiles, function(x) {
    x = tools::file_path_sans_ext(basename(x))
    if (endsWith(x, ".tmad")) {
        x = substr(x, 1, nchar(x) - 5)
    }
    if (!is.null(sample_name_fun)) {
        x = sample_name_fun(x)
    }
    x
})

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

data = data.frame(Sample = sams, tMAD = tmads)
if (is.character(metafile) && file.exists(metafile) && length(meta_cols) > 0) {
    metadf = read.table(metafile, header=T, row.names=NULL, sep="\t", stringsAsFactors=F)
    if (!is.null(metadf$Sample)) {
        metadf$Sample = as.character(metadf$Sample)
    } else {
        colnames(metadf)[1] = "Sample"
    }
    meta = metadf[, c("Sample", meta_cols), drop=FALSE]
    colnames(meta) = c("Sample", meta_cols)
    data = data %>% left_join(meta, by="Sample")
}

# save the data
write.table(data, file=file.path(outdir, "tMAD.txt"), sep="\t", quote=F, row.names=F)

# bar plot for all samples without grouping
p <- BarPlot(
    data = data,
    x = "Sample",
    y = "tMAD",
    x_text_angle = 90
)

png(file.path(outdir, "tMAD.png"), width=400 + length(sams) * 12, height=800, res=100)
print(p)
dev.off()

# Do it for each group col
if (!is.null(group_cols)) {
    for (group_col in group_cols) {
        if (!grepl(",", group_col, fixed = TRUE)) {
            # Bar plot with this group_col, but with different fill colors
            # for each group, and samples from the same group are next to each other
            gdata <- data %>% arrange(!!sym(group_col)) %>% mutate(Sample=factor(Sample, levels=unique(Sample)))
            p <- BarPlot(
                data = gdata,
                x = "Sample",
                y = "tMAD",
                fill = group_col,
                x_text_angle = 90
            )

            png(file.path(outdir, paste0("tMAD_", group_col, "_bar.png")), width=400 + length(sams) * 12, height=600, res=100)
            print(p)
            dev.off()

            # Box plot overlays with violin plot with this group_col
            p <- ViolinPlot(
                data = gdata,
                x = group_col,
                y = "tMAD",
                x_text_angle = 90,
                add_box = TRUE,
                add_point = TRUE,
                comparisons = TRUE,
                sig_label = "p.format"
            )

            png(file.path(outdir, paste0("tMAD_", group_col, "_box_violin.png")), width=1000, height=600, res=100)
            print(p)
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
            gdata <- data %>% arrange(!!sym(group_col1), !!sym(group_col2)) %>%
                mutate(Sample=factor(Sample, levels=unique(Sample)))
            p <- BarPlot(
                data = gdata,
                x = "Sample",
                y = "tMAD",
                split_by = group_col1,
                fill = group_col2,
                x_text_angle = 90,
                ncol = 2
            )

            png(
                file.path(outdir, paste0("tMAD_", group_col, "_bar.png")),
                width=1000,
                height=length(unique(data[[group_col1]])) * 200,
                res=100
            )
            print(p)
            dev.off()

            # Do the same for Voilin + boxplot
            p <- ViolinPlot(
                data = gdata,
                x = group_col2,
                y = "tMAD",
                split_by = group_col1,
                x_text_angle = 90,
                add_box = TRUE,
                add_point = TRUE,
                comparisons = TRUE,
                sig_label = "p.format",
                ncol = 2
            )

            png(
                file.path(outdir, paste0("tMAD_", group_col, "_box_violin.png")),
                width=1000,
                height=length(unique(data[[group_col1]])) * 200,
                res=100
            )
            print(p)
            dev.off()
        }
    }
}
