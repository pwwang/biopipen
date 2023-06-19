library(ggplot2)
library(ggprism)
library(dplyr)
library(tidyr)
library(tibble)

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

data = data.frame(Sample = sams, tMAD = tmads)
if (file.exists(metafile) && length(group_cols) > 0) {
    metadf = read.table(metafile, header=T, row.names=NULL, sep="\t", stringsAsFactors=F)
    sample_col = colnames(metadf)[1]
    meta = metadf[, c(sample_col, group_cols), drop=FALSE]
    colnames(meta) = c("Sample", group_cols)
    data = data %>% left_join(meta, by="Sample")
}

# save the data
write.table(data, file=file.path(outdir, "tMAD.txt"), sep="\t", quote=F, row.names=F)

# bar plot for all samples without grouping
p = ggplot(data, aes(x=Sample, y=tMAD)) +
    geom_bar(stat="identity", fill="steelblue") +
    theme_prism() +
    theme(
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position = "none",
    ) +
    labs(
        x = NULL,
        y = "tMAD",
    )

png(file.path(outdir, "tMAD.png"), width=400 + length(sams) * 10, height=800, res=100)
print(p)
dev.off()

# Do it for each group col
if (!is.null(group_cols)) {
    for (group_col in group_cols) {
        # Bar plot with this group_col, but with different fill colors
        # for each group, and samples from the same group are next to each other
        p = ggplot(
                data %>% arrange(!!sym(group_col)) %>% mutate(Sample=factor(Sample, levels=Sample)),
                aes(x=Sample, y=tMAD, fill=!!sym(group_col))
            ) +
            geom_bar(stat="identity") +
            theme_prism() +
            theme(
                axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size=12),
                axis.text.y = element_text(size=12),
            ) +
            labs(
                x = NULL,
                y = "tMAD",
            )

        png(file.path(outdir, paste0("tMAD_", group_col, "_bar.png")), width=400 + length(sams) * 10, height=600, res=100)
        print(p)
        dev.off()

        # Box plot overlays with violin plot with this group_col
        p = ggplot(data, aes(x=!!sym(group_col), y=tMAD)) +
            geom_boxplot(outlier.shape=NA, fill="white", color="black") +
            geom_violin(fill="steelblue", alpha=0.5) +
            theme_prism() +
            theme(
                axis.title.x = element_text(size=12),
                axis.title.y = element_text(size=12),
                axis.text.y = element_text(size=12),
            ) +
            labs(
                x = group_col,
                y = "tMAD",
            )

        png(file.path(outdir, paste0("tMAD_", group_col, "_box_violin.png")), width=1000, height=600, res=100)
        print(p)
        dev.off()
    }
}
