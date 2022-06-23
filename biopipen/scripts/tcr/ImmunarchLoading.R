source("{{biopipen_dir}}/utils/misc.R")

# Loading 10x data into immunarch
library(immunarch)
library(dplyr)
library(tidyr)
library(tibble)
library(glue)
library(bracer)

metafile = {{ in.metafile | quote }}
rdsfile = {{ out.rdsfile | quote }}
metatxt = {{ out.metatxt | quote }}
tmpdir = {{ envs.tmpdir | quote }}
mode = {{ envs.mode | quote }}
metacols = {{ envs.metacols | r}}

metadata = read.table(
    metafile,
    header = TRUE,
    row.names = NULL,
    sep = "\t",
    check.names = FALSE
)

if (!"Sample" %in% colnames(metadata)) {
    stop("Error: Column `Sample` is not found in metafile.")
}
if (!"TCRDir" %in% colnames(metadata)) {
    stop("Error: Column `TCRDir` is not found in metafile.")
}

## --------------------------------------------------
## Helpers

get_contig_annofile = function(dir, sample, warn=TRUE) {
    annofilepat = paste0(
        "*", "{all,filtered}", "_contig_annotations.csv*"  # .gz
    )
    annofiles = glob(file.path(as.character(dir), annofilepat))
    if (length(annofiles) == 0) {
        stop(paste(
            "Cannot find neither `filtered_contig_annotations.csv[.gz]` nor",
            "`all_contig_annotations.csv[.gz]`",
            "in given TCRDir for sample:",
            sample
        ))
    } else if (length(annofiles) > 1) {
        if (warn) {
            warning(paste(
                "Found more than one file in given TCRDir for sample:",
                sample
            ))
        }
        for (annofile in annofiles) {
            # use filtered if both filtered_ and all_ are found
            if (grepl("filtered", annofile)) {
                annofiles = annofile
                break
            }
            # give a warning if only all_ is found
            if (warn) {
                warning(paste(
                    "Using all_contig_annotations as",
                    "filtred_config_annotations not found",
                    "in given TCRDir for sample:",
                    sample
                ))
            }
        }
    }
    annofiles[1]
}

datadir = tempfile(pattern = "immunarch-", tmpdir = tmpdir)
dir.create(datadir, showWarnings = FALSE)

# Find filtered_contig_annotations.csv and link then in datadir
for (i in seq_len(nrow(metadata))) {
    sample = as.character(metadata[i, "Sample"])

    annofile = get_contig_annofile(metadata[i, "TCRDir"], sample)
    filename = basename(annofile)
    stem = sub("\\.gz$", "", filename)
    stem = sub("\\.csv", "", stem)
    ext = substr(filename, nchar(stem) + 1, nchar(filename))
    file.symlink(normalizePath(annofile), file.path(datadir, paste0(sample, ext)))
}

immdata = repLoad(datadir, .mode=mode)
if (mode == "single") {
    data = immdata$data
    immdata$tra = list()
    immdata$multi = list()
    immdata$data = list()
    for (name in names(data)) {
        if (endsWith(name, "_TRA")) {
            immdata$tra[substr(name, 1, nchar(name) - 4)] = data[name]
        } else if (endsWith(name, "_TRB")) {
            immdata$data[substr(name, 1, nchar(name) - 4)] = data[name]
        } else {  # "X_Multi"
            immdata$multi[substr(name, 1, nchar(name) - 6)] = data[name]
        }

    }
}

if (mode == "single") {
    immdata$meta  = immdata$meta |>
        filter(endsWith(Sample, "_TRB")) |>
        mutate(Sample = substr(Sample, 1, nchar(Sample) - 4)) |>
        select(-"Source")
}
immdata$meta = left_join(
    immdata$meta,
    metadata %>% select(-"TCRDir"),
    by = "Sample"
)

saveRDS(immdata, file=rdsfile)

metadf = do_call(rbind, lapply(seq_len(nrow(immdata$meta)), function(i) {
    # Clones  Proportion   CDR3.aa                       Barcode
    # 5      4 0.008583691 CAVRDTGNTPLVF;CASSEYSNQPQHF   GTTCGGGCACTTACGA-1;TCTCTAAGTACCAGTT-1
    # 6      4 0.008583691 CALTQAAGNKLTF;CASRPEDLRGQPQHF GCTTGAAGTCGGCACT-1;TACTCGCTCCTAAGTG-1
    cldata = immdata$data[[i]][, unique(c(metacols, "Barcode"))]
    # # A tibble: 4 × 5
    # Sample                  Patient     Timepoint Tissue
    # <chr>                   <chr>       <chr>     <chr>
    # 1 MC1685Pt011-Baseline-PB MC1685Pt011 Baseline  PB
    mdata = as.list(immdata$meta[i, ])
    for (mname in names(mdata)) {
        assign(mname, mdata[[mname]])
    }

    cldata |>
        separate_rows(Barcode, sep=";") |>
        distinct(Barcode, .keep_all = TRUE) |>
        mutate(Barcode = glue("{{envs.prefix}}{Barcode}")) |>
        column_to_rownames("Barcode")

}))
write.table(metadf, metatxt, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
