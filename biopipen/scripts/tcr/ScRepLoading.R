library(rlang)
library(bracer)
library(scRepertoire)
library(biopipen.utils)

metafile <- {{in.metafile | r}}
outfile <- {{out.outfile | r}}
combineTCR_args <- {{envs.combineTCR | r}}
combineBCR_args <- {{envs.combineBCR | r}}
type <- {{envs.type | r}}
exclude <- {{envs.exclude | r}}
format <- {{envs.format | r}}
tmpdir <- {{envs.tmpdir | r}}

type = toupper(type)
if (length(exclude) == 1) {
    exclude <- strsplit(exclude, ",")[[1]]
}

log <- get_logger()

log$info("Loading metadata ...")
metadata <- read.table(metafile, header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE)
if (type == "AUTO") {
    if ("TCRData" %in% colnames(metadata) && "BCRData" %in% colnames(metadata)) {
        log$warn("Both TCRData and BCRData columns found in metadata. Defaulting to TCR.")
        type <- "TCR"
    } else if ("TCRData" %in% colnames(metadata)) {
        type <- "TCR"
    } else if ("BCRData" %in% colnames(metadata)) {
        type <- "BCR"
    } else {
        stop("Error: Neither TCRData nor BCRData column found in metadata.")
    }
}

data_column <- ifelse(type == "TCR", "TCRData", "BCRData")
combine_fn <- ifelse(type == "TCR", combineTCR, combineBCR)
combine_args <- if (type == "TCR") { combineTCR_args } else { combineBCR_args }

stopifnot("Error: Column `Sample` is not found in metafile." = "Sample" %in% colnames(metadata))
if (!data_column %in% colnames(metadata)) {
    stop(paste0("Error: Column `", data_column, "` is not found in metafile."))
}
rownames(metadata) <- metadata$Sample

.gunzip <- function(input, output) {
    # Open connections
    con_in <- gzfile(input, "rt")   # "rt" = read text mode
    con_out <- file(output, "wt")  # "wt" = write text mode

    # Read line by line and write
    while(length(line <- readLines(con_in, n = 10, warn = FALSE)) > 0) {
        writeLines(line, con_out)
    }

    # Close connections
    close(con_in)
    close(con_out)
}

get_file_name <- function(fmt) {
    if (is.null(fmt)) { return("filtered_contig_annotations.csv") }
    fmt <- tolower(fmt)
    if (fmt == "10x") { return("filtered_contig_annotations.csv") }
    if (fmt == "airr") { return("airr_rearrangement.tsv") }
    if (fmt == "bd") { return("Contigs_AIRR.tsv") }
    if (fmt == "dandelion") { return("all_contig_dandelion.tsv") }
    if (fmt == "immcantation") { return("data.tsv") }
    if (fmt == "json") { return("contigs.json") }
    if (fmt == "parsebio") { return("barcode_report.tsv") }
    if (fmt == "mixcr") { return("clones.tsv") }
    if (fmt == "omniscope") { return("contigs.csv") }
    if (fmt == "trust4") { return("barcode_report.tsv") }
    if (fmt == "wat3r") { return("barcode_results.csv") }

    stop("Unsupported format: ", fmt)
}

get_format <- function(filename) {
    if (identical(filename, "filtered_contig_annotations.csv")) { return("10X") }
    if (identical(filename, "airr_rearrangement.tsv")) { return("AIRR") }
    if (identical(filename, "Contigs_AIRR.tsv")) { return("BD") }
    if (identical(filename, "all_contig_dandelion.tsv")) { return("Dandelion") }
    if (identical(filename, "data.tsv")) { return("Immcantation") }
    if (endsWith(filename, ".json")) { return("JSON") }
    if (identical(filename, "barcode_report.tsv")) { return("ParseBio") }
    if (identical(filename, "clones.tsv")) { return("MiXCR") }
    if (identical(filename, "contigs.csv")) { return("Omniscope") }
    # if (identical(filename, "barcode_report.tsv")) { return("TRUST4") }
    if (identical(filename, "barcode_results.csv")) { return("WAT3R") }

    return("10X")
}

# helper function
get_contig_dir <- function(input, sample, fmt) {
    if (is.na(input) || !is.character(input) || nchar(input) == 0 || input == "NA") {
        warning(paste0("No path found for sample: ", sample), immediate. = TRUE)
        return(list(NULL, fmt))
    }
    if (!file.exists(input)) {
        stop(paste0("Input path does not exist for sample: ", sample, ": ", input))
    }
    if (dir.exists(input)) {
        return(list(input, fmt))
    }
    # file
    filedir <- file.path(tmpdir, slugify(sample))
    dir.create(filedir, recursive = TRUE, showWarnings = FALSE)

    # if it is gzipped
    if (grepl("\\.gz$", input)) {
        flatfile <- file.path(filedir, sub("\\.gz$", "", basename(input)))
        .gunzip(input, flatfile)
        input <- flatfile
    }

    fmt <- fmt %||% get_format(basename(input))
    filename <- get_file_name(fmt)
    file.symlink(input, file.path(filedir, filename))

    return(list(filedir, fmt))
}

load_contig <- function(input, sample, fmt) {
    log$info("- Sample: {sample}")
    dirfmt <- get_contig_dir(input, sample, fmt)
    dir <- dirfmt[[1]]
    fmt <- dirfmt[[2]]
    if (is.null(dir)) { return(NULL) }
    x <- loadContigs(dir, format = fmt %||% "10X")
    x <- x[[1]]
    x$sample <- NULL
    if (identical(fmt %||% "10X", "10X") && colnames(x)[1] == "X") {
        x$X <- NULL
    }

    x
}


log$info("Reading {type} data ...")
contig_list <- lapply(seq_len(nrow(metadata)), function(i) {
    sample <- as.character(metadata$Sample[i])
    path <- metadata[[data_column]][i]
    load_contig(path, sample, fmt = format)
})
names(contig_list) <- as.character(metadata$Sample)
contig_list <- contig_list[!sapply(contig_list, is.null)]

log$info("Combining {type} data and adding meta data ...")
if (isTRUE(combine_args$samples)) {
    combine_args$samples <- names(contig_list)
}
combine_args$input.data <- contig_list
screp_data <- do_call(combine_fn, combine_args)
for (col in colnames(metadata)) {
    if (col %in% exclude) { next }
    screp_data <- addVariable(screp_data, col, metadata[names(screp_data), col])
}

rm(contig_list, combine_args)

log$info("Saving {type} data ...")
save_obj(screp_data, outfile)
