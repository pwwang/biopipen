library(parallel)
library(dplyr)
library(biopipen.utils)
library(CNAclinic)

# https://github.com/sdchandra/CNAclinic/issues/4
.reorderByChrom.patched <- function(x){
    chromosome <- as.character(x$chromosome)
    chromosome[which(chromosome == "X")] <- "23"
    chromosome[which(chromosome == "Y")] <- "24"
    chromosome[which(chromosome == "MT")] <- "25"

    x$chromosome <- as.numeric(chromosome)
    # Error in xtfrm.data.frame(x) : cannot xtfrm data frames
    # x <- x[order(x["chromosome"], x["start"]), ]
    x <- x[order(x[, "chromosome"], x[, "start"]), ]

    x$chromosome <- as.character(x$chromosome)
    # Replace 23 by X:
    x$chromosome[which(x$chromosome == "23")] <- "X"

    # Replace 24 by Y
    x$chromosome[which(x$chromosome == "24")] <- "Y"

    # Replace 25 by MT
    x$chromosome[which(x$chromosome == "25")] <- "MT"

    return(x)
}

monkey_patch("CNAclinic", ".reorderByChrom", .reorderByChrom.patched)

metafile = {{in.metafile | r}}
outdir = {{out.outdir | r}}
ncores = {{envs.ncores | int}}
binsizer = {{envs.binsizer | r}}
binsize = {{envs.binsize | r}}
seed = {{envs.seed | int}}
genome = {{envs.genome | r}}
run_args = {{envs.run_args | r}}
plot_args = {{envs.plot_args | r}}
plot_multi_args = {{envs.plot_multi_args | r}}

set.seed(seed)

metadata = read.table(metafile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metacols = colnames(metadata)

# check if the metadata file has the required columns
if (!("Bam" %in% metacols)) {
    stop("The metadata file must have a column named 'Bam'")
}

if (("Group" %in% metacols) && !("Patient" %in% metacols)) {
    stop("The metadata file must have a column named 'Patient' if it has a column named 'Group'")
}

if (!("Binsizer" %in% metacols) && is.null(binsizer) && is.null(binsize)) {
    stop(
        "The metadata file must have a column named 'Binsizer' or ",
        "the `envs.binsizer` must be specified when no `envs.binsize` is provided. ",
        "The Binsizer column should indicate which samples are to be used for binsize selection."
    )
}

# add missing columns
if (!("Sample" %in% metacols)) {
    metadata$Sample = tools::file_path_sans_ext(basename(metadata$Bam))
}

if (!("Binsizer" %in% metacols) && is.null(binsize)) {
    if (!is.atomic(binsizer)) {
        metadata$Binsizer = metadata$Sample %in% binsizer
    } else if (!is.numeric(binsizer)) {
        stop("The `envs.binsizer` must be a numeric value or a vector of sample names")
    } else {
        if (binsizer > 1) {binsizer = binsizer / nrow(metadata)}

        metadata$Binsizer = sample(
            c(TRUE, FALSE),
            nrow(metadata),
            replace = TRUE,
            prob = c(binsizer, 1 - binsizer)
        )
    }
}

# convert Binsizer to logical
if (!is.logical(metadata$Binsizer) && is.null(binsize)) {
    metadata$Binsizer = metadata$Binsizer %in% c("TRUE", "T", "1", "Yes", "Y", "True")
}

if (is.null(binsize)) {
    binsizer_df = metadata[metadata$Binsizer, ]
    binsizer_bam = binsizer_df$Bam
    binsizer_sample = binsizer_df$Sample

    plotList = optimalBinsize(binsizer_bam, binsizer_sample)

    binsize_plot = file.path(outdir, "binsize_plot.png")
    png(binsize_plot, width = 1000, height = 1000, res = 100)
    print(plotList[[1]])
    dev.off()

    y = plotList[[1]]$layers[[2]]$data$yintercept
    bsdata = plotList[[1]]$data |> filter(between(value, y[1], y[2]))
    binsize = bsdata |>
        count(binsize) |>
        arrange(desc(n)) |>
        head(1) |>
        pull(binsize) |>
        as.character() |>
        as.numeric()

    write.table(bsdata, file = file.path(outdir, "binsizes.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}
cat("The optimal binsize is ", binsize, "Kbp", file = file.path(outdir, "selected_binsize.txt"))

# Process BAM files
bamfiles = metadata$Bam
samples = metadata$Sample

do_one_sample = function(i) {
    bamfile = bamfiles[i]
    sample = samples[i]
    refSamples = NULL

    if ("Group" %in% metacols) {
        group = metadata |> filter(Sample == sample) |> pull(Group)
        if (group == "Control") {
            return (NULL)
        } else {
            patient = metadata |> filter(Sample == sample) |> pull(Patient)
            refSamples = metadata |> filter(Patient == patient | Group == "Control") |> pull(Sample)
        }
    }

    processedData = processForSegmentation(
        bamfile,
        sample,
        refSamples=refSamples,
        binSize=binsize / 1000
    )

    run_args_i = run_args
    run_args_i$x = processedData
    run_args_i$genome = genome
    CNAData = do_call(runSegmentation, run_args_i)

    plot_args_i = plot_args
    plot_args_i$object = CNAData
    genomewide_plot <- tryCatch({
        do_call(plotSampleData, plot_args_i)
    }, error = function(e) {
        message("Error in plotting genomewide data for sample ", sample, ": ", e$message)
        return(ggplot2::ggplot() + ggplot2::labs(title = paste("Error in plotting genomewide data for sample", sample)))
    })

    odir = file.path(outdir, sample)
    dir.create(odir, recursive = TRUE, showWarnings = FALSE)
    png(file.path(odir, "genomewide.png"), width = 1000, height = 250, res = 100)
    print(genomewide_plot[[1]])
    dev.off()

    # save the CNAData
    exportData(CNAData,
               dataType="calls", fileType="igv",
               fileName=file.path(odir, "calls.igv"))

    # for a list of genes as a csv file
    # exportData(combined,
    #         geneInfo=geneInfo,
    #         dataType="summary", fileType="csv",
    #         fileName=file.path(odir, "Genes_segSummary.csv"))
    if (is.list(plot_multi_args)) {
        return(CNAData)
    } else {
        return(NULL)
    }
}

if (ncores > 1) {
    res = mclapply(1:length(bamfiles), do_one_sample, mc.cores = ncores)
} else {
    res = lapply(1:length(bamfiles), do_one_sample)
}

if (is.list(plot_multi_args)) {
    combined = NULL
    for (i in 1:length(res)) {
        if (!is.null(combined)) {
            combined = combineData(combined, res[[i]])
        } else {
            combined = res[[i]]
        }
    }

    plot_multi_args$object = combined
    multiplot = do_call(plotMultiSampleData, plot_multi_args)
    png(file.path(outdir, "multiplot.png"), width = 1000, height = 1000, res = 100)
    print(multiplot)
    dev.off()
}
