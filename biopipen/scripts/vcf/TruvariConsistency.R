{{ biopipen_dir | joinpaths: "utils", "plot.R" | source_r }}
# to compile the expressions
library(ComplexHeatmap)

vcfs = {{in.vcfs | r}}
outdir = {{out.outdir | r}}
truvari = {{envs.truvari | r}}
heatmap = {{envs.heatmap | r}}

consistency_file = file.path(outdir, "consistency.txt")
consistency_plot = file.path(outdir, "consistency.png")

system2(
    paste(c(truvari, "consistency", vcfs)),
    stdout = consistency_file
)

if (isFALSE(heatmap)) {
    quit(save = "no")
}

# parse the consistency file to get
# 1. The sample name
# 2. The presence matrix

read_consistency = function() {
    samples = c()
    presence = list()
    con = file(consistency_file, "r")

    SAMPLES_BLOCK = FALSE
    PRESENCE_BLOCK = FALSE
    CNV_INDEX = 1

    while (TRUE) {
        line = readLines(con, 1)
        if (length(line) == 0) break

        if (line == "#File	NumCalls") {
            SAMPLES_BLOCK = TRUE
        } else if (line == "#Group	Total	TotalPct	PctOfFileCalls") {
            PRESENCE_BLOCK = TRUE
        } else if (line == "#") {
            SAMPLES_BLOCK = FALSE
            PRESENCE_BLOCK = FALSE
        } else if (SAMPLES_BLOCK) {
            samples = c(
                samples,
                tools::file_path_sans_ext(basename(strsplit(line, "\t")[[1]][1]))
            )
        } else if (PRESENCE_BLOCK) {
            presence[[paste0("CNV_", CNV_INDEX)]] = as.integer(strsplit(
                strsplit(line, "\t")[[1]][1], split=""
            )[[1]])
            CNV_INDEX = CNV_INDEX + 1
        }
    }
    close(con)

    out = as.data.frame(presence)
    rownames(out) = samples
    #          CNV_1   CNV_2   CNV_3   ...
    # sample1  1       0       0       ...
    # sample2  0       1       0       ...
    # sample3  0       0       1       ...
    # ...
    t(out)
}

df = read_consistency()

if (is.null(heatmap$args)) {
    heatmap$args = list()
}

if (is.null(heatmap$draw)) {
    heatmap$draw = list()
}

if (!is.null(heatmap$annofile)) {
    annos = read.table(
        heatmap$annofile,
        header = T,
        row.names = 1,
        sep = "\t"
    )
    ha = HeatmapAnnotation(df = annos)
    heatmap$args$top_annotation = ha
}

plotHeatmap(
    df,
    args = heatmap$args,
    draw = heatmap$draw,
    devpars = heatmap$devpars,
    outfile = consistency_plot
)
