{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}

library(rlang)
library(bracer)
library(scRepertoire)

metafile <- {{in.metafile | quote}}
outfile <- {{out.outfile | quote}}
combineTCR_args <- {{envs.combineTCR | r}}
exclude <- {{envs.exclude | r}}
if (length(exclude) == 1) {
    exclude <- strsplit(exclude, ",")[[1]]
}

log_info("Loading metadata ...")
metadata <- read.table(metafile, header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE)
stopifnot("Error: Column `Sample` is not found in metafile." = "Sample" %in% colnames(metadata))
stopifnot("Error: Column `TCRData` is not found in metafile." = "TCRData" %in% colnames(metadata))
rownames(metadata) <- metadata$Sample

# helper function
get_contig_annofile <- function(dir, sample, warn = TRUE) {
    if (is.na(dir) || !is.character(dir) || nchar(dir) == 0 || dir == "NA") {
        warning(paste0("No path found for sample: ", sample), immediate. = TRUE)
        return (NULL)
    }
    if (file.exists(dir) && !dir.exists(dir)) {
        return(dir)
    }

    annofilepat <- paste0("*", "{all,filtered}", "_contig_annotations.csv*")  # .gz
    annofiles <- glob(file.path(as.character(dir), annofilepat))
    if (length(annofiles) == 0) {
        stop(
            "Cannot find neither `filtered_contig_annotations.csv[.gz]` nor",
            "`all_contig_annotations.csv[.gz]` in given TCRData for sample: ",
            sample
        )
    }
    if (length(annofiles) > 1) {
        if (warn) {
            warning("Found more than one file in given TCRData for sample: ", sample, immediate. = TRUE)
        }
        for (annofile in annofiles) {
            # use filtered if both filtered_ and all_ are found
            if (grepl("filtered", annofile)) {
                annofiles <- annofile
                break
            }
            # give a warning if only all_ is found
            if (warn) {
                warning("Using all_contig_annotations as filtred_config_annotations not found ",
                        "in given TCRData for sample: ", sample,
                        immediate. = TRUE
                )
            }
        }
    }
    annofiles[1]
}

# for (i in seq_len(nrow(metadata))) {
#     sample <- as.character(metadata$Sample[i])
#     annofile <- get_contig_annofile(metadata$TCRData[i], sample)
#     if (is.null(annofile)) { next }

#     anno <- read.delim2(annofile, sep = ",", header = TRUE, stringsAsFactors = FALSE)
#     # Add cdr1, cdr2, fwr1, fwr2, etc columns
#     anno$cdr1 <- anno$cdr1 %||% ""
#     anno$cdr1_nt <- anno$cdr1_nt %||% ""
#     anno$cdr2 <- anno$cdr2 %||% ""
#     anno$cdr2_nt <- anno$cdr2_nt %||% ""
#     anno$fwr1 <- anno$fwr1 %||% ""
#     anno$fwr1_nt <- anno$fwr1_nt %||% ""
#     anno$fwr2 <- anno$fwr2 %||% ""
#     anno$fwr2_nt <- anno$fwr2_nt %||% ""
#     anno$fwr3 <- anno$fwr3 %||% ""
#     anno$fwr3_nt <- anno$fwr3_nt %||% ""
#     anno$fwr4 <- anno$fwr4 %||% ""
#     anno$fwr4_nt <- anno$fwr4_nt %||% ""

#     annotfile = file.path(datadir, paste0(sample, ".csv"))
#     write.table(anno, annotfile, sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
# }

log_info("Reading TCR data ...")
contig_list <- lapply(seq_len(nrow(metadata)), function(i) {
    sample <- as.character(metadata$Sample[i])
    annofile <- get_contig_annofile(metadata$TCRData[i], sample)
    if (is.null(annofile)) { return (NULL) }

    log_info("- Sample: {sample} ...")
    anno <- read.delim2(annofile, sep = ",", header = TRUE, stringsAsFactors = FALSE)
    # Add cdr1, cdr2, fwr1, fwr2, etc columns for compatibility
    anno$cdr1 <- anno$cdr1 %||% ""
    anno$cdr1_nt <- anno$cdr1_nt %||% ""
    anno$cdr2 <- anno$cdr2 %||% ""
    anno$cdr2_nt <- anno$cdr2_nt %||% ""
    anno$fwr1 <- anno$fwr1 %||% ""
    anno$fwr1_nt <- anno$fwr1_nt %||% ""
    anno$fwr2 <- anno$fwr2 %||% ""
    anno$fwr2_nt <- anno$fwr2_nt %||% ""
    anno$fwr3 <- anno$fwr3 %||% ""
    anno$fwr3_nt <- anno$fwr3_nt %||% ""
    anno$fwr4 <- anno$fwr4 %||% ""
    anno$fwr4_nt <- anno$fwr4_nt %||% ""

    anno
})
names(contig_list) <- as.character(metadata$Sample)
contig_list <- contig_list[!sapply(contig_list, is.null)]

log_info("Combining TCR data and adding meta data ...")
if (isTRUE(combineTCR_args$samples)) {
    combineTCR_args$samples <- names(contig_list)
}
combineTCR_args$input.data <- contig_list
screp_data <- do_call(combineTCR, combineTCR_args)
for (col in colnames(metadata)) {
    if (col %in% exclude) { next }
    screp_data <- addVariable(screp_data, col, metadata[names(screp_data), col])
}

rm(contig_list, combineTCR_args)

log_info("Saving TCR data ...")
saveRDS(screp_data, outfile)
