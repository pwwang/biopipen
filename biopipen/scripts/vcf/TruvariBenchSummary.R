library(rjson)
library(rlang)
library(dplyr)
library(plotthis)

indirs = {{in.indirs | r}}
outdir = {{out.outdir | r}}
plots = {{envs.plots | r}}
devpars = {{envs.devpars | r}}

read_summary = function() {
    # Read the summary file from each directory

    summaries = NULL
    for (indir in indirs) {
        summary = fromJSON(file=file.path(indir, "summary.json"))
        summary$gt_matrix = NULL
        summary$Sample = sub(".truvari_bench", "", basename(indir), fixed=T)
        summaries = bind_rows(summaries, summary)
    }
    summaries
}

get_devpars = function() {
    if (!is.null(devpars)) {
        return (devpars)
    }
    n_samples = length(indirs)
    list(
        res = 100,
        height = 1000,
        width = 100 * n_samples + 150
    )
}

plot_summary = function(col) {
    outfile = file.path(outdir, paste0(col, ".png"))
    p <- plotthis::BarPlot(
        summaries,
        x = "Sample",
        y = col,
        x_text_angle = 90
    )
    devpars <- get_devpars()
    png(
        filename = outfile,
        width = devpars$width,
        height = devpars$height,
        res = devpars$res
    )
    print(p)
    dev.off()
}

main = function() {
    for (col in plots) {
        plot_summary(col)
    }
}

{
    summaries = read_summary()
    write.table(
        summaries %>% select(Sample, everything()),
        file=file.path(outdir, "summary.txt"),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )
    main()
}
