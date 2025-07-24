{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "single_cell.R" | source_r }}

# Loading 10x data into immunarch
library(immunarch)
library(dplyr)
library(tidyr)
library(tibble)
library(glue)
library(bracer)

metafile = {{ in.metafile | r }}
rdsfile = {{ out.rdsfile | r }}
metatxt = {{ out.metatxt | r }}
tmpdir = {{ envs.tmpdir | r }}
mode = {{ envs.mode | r }}
extracols = {{ envs.extracols | r}}
prefix = {{ envs.prefix | r }}

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
if (!"TCRData" %in% colnames(metadata)) {
    stop("Error: Column `TCRData` is not found in metafile.")
}

## --------------------------------------------------
## Helpers

get_contig_annofile = function(dir, sample, warn=TRUE) {
    if (is.na(dir) || !is.character(dir) || nchar(dir) == 0 || dir == "NA") {
        warning(paste0("No path found for sample: ", sample), immediate. = TRUE)
        return (NULL)
    }

    annofilepat = paste0(
        "*", "{all,filtered}", "_contig_annotations.csv*"  # .gz
    )
    annofiles = glob(file.path(as.character(dir), annofilepat))
    if (length(annofiles) == 0) {
        stop(paste(
            "Cannot find neither `filtered_contig_annotations.csv[.gz]` nor",
            "`all_contig_annotations.csv[.gz]`",
            "in given TCRData for sample:",
            sample
        ))
    } else if (length(annofiles) > 1) {
        if (warn) {
            warning(
                paste(
                    "Found more than one file in given TCRData for sample:",
                    sample
                ),
                immediate. = TRUE
            )
        }
        for (annofile in annofiles) {
            # use filtered if both filtered_ and all_ are found
            if (grepl("filtered", annofile)) {
                annofiles = annofile
                break
            }
            # give a warning if only all_ is found
            if (warn) {
                warning(
                    paste(
                        "Using all_contig_annotations as",
                        "filtred_config_annotations not found",
                        "in given TCRData for sample:",
                        sample
                    ),
                    immediate. = TRUE
                )
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

    annofile = get_contig_annofile(metadata[i, "TCRData"], sample)
    if (is.null(annofile)) {
        next
    }
    anno = read.delim2(annofile, sep=",", header=TRUE, stringsAsFactors=FALSE)
    # Add cdr1, cdr2, fwr1, fwr2, etc columns
    if (!"cdr1" %in% colnames(anno)) {
        anno$cdr1 = ""
    }
    if (!"cdr1_nt" %in% colnames(anno)) {
        anno$cdr1_nt = ""
    }
    if (!"cdr2" %in% colnames(anno)) {
        anno$cdr2 = ""
    }
    if (!"cdr2_nt" %in% colnames(anno)) {
        anno$cdr2_nt = ""
    }
    if (!"fwr1" %in% colnames(anno)) {
        anno$fwr1 = ""
    }
    if (!"fwr1_nt" %in% colnames(anno)) {
        anno$fwr1_nt = ""
    }
    if (!"fwr2" %in% colnames(anno)) {
        anno$fwr2 = ""
    }
    if (!"fwr2_nt" %in% colnames(anno)) {
        anno$fwr2_nt = ""
    }
    if (!"fwr3" %in% colnames(anno)) {
        anno$fwr3 = ""
    }
    if (!"fwr3_nt" %in% colnames(anno)) {
        anno$fwr3_nt = ""
    }
    if (!"fwr4" %in% colnames(anno)) {
        anno$fwr4 = ""
    }
    if (!"fwr4_nt" %in% colnames(anno)) {
        anno$fwr4_nt = ""
    }
    annotfile = file.path(datadir, paste0(sample, ".csv"))
    write.table(anno, annotfile, sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

    # filename = basename(annofile)
    # stem = sub("\\.gz$", "", filename)
    # stem = sub("\\.csv", "", stem)
    # ext = substr(filename, nchar(stem) + 1, nchar(filename))
    # file.symlink(normalizePath(annofile), file.path(datadir, paste0(sample, ext)))
}

log_info("Loading TCR data ...")
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
    immdata$meta  = immdata$meta %>%
        filter(endsWith(Sample, "_TRB")) %>%
        mutate(Sample = substr(Sample, 1, nchar(Sample) - 4)) %>%
        select(-"Source")
}
immdata$meta = left_join(
    immdata$meta,
    metadata %>% select(-"TCRData"),
    by = "Sample"
)

immdata$prefix = prefix

saveRDS(immdata, file=rdsfile)

log_info("Saving cell-level data ...")
exdata <- expand_immdata(immdata, cell_id = "Barcode") %>%
    distinct(Sample, Barcode, .keep_all = TRUE) %>%
    mutate(Barcode = glue(paste0(prefix, "{Barcode}"))) %>%
    select(any_of(c(
        colnames(immdata$meta),
        "Barcode",
        "CDR3.aa",
        "Clones",
        "Proportion",
        "V.name",
        "D.name",
        "J.name",
        extracols
    ))) %>%
    column_to_rownames("Barcode")

write.table(exdata, metatxt, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
