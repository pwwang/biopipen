# Loading 10x data into immunarch
library(immunarch)
library(dplyr)
library(tidyr)
library(bracer)

metafile = {{ in.metafile | quote }}
rdsfile = {{ out.rdsfile | quote }}
tmpdir = {{ envs.tmpdir | quote }}

metadata = read.table(
    metafile,
    header = TRUE,
    row.names = NULL,
    sep = "\t",
    check.names = FALSE
)

meta_cols = colnames(metadata)
if (!"Sample" %in% meta_cols) {
    stop("Error: Column `Sample` is not found in metafile.")
}
if (!"TCRDir" %in% meta_cols) {
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

immdata = repLoad(datadir)
# drop TRAs for paired data
immdata$single = list()
immdata$raw = list()
for (sample in names(immdata$data)) {
    annofile = get_contig_annofile(
        metadata[metadata$Sample == sample, "TCRDir"],
        sample,
        warn=FALSE
    )
    immdata$raw[[sample]] = read.csv2(
        annofile,
        header = TRUE,
        row.names = NULL,
        sep = ",",
        check.names = FALSE
    )
    immdata$single[[sample]] = immdata$data[[sample]] %>%
        separate_rows(
            CDR3.nt, CDR3.aa, V.name, D.name, J.name, Sequence, chain,
            sep=";"
        ) %>%
        filter(chain != "TRA") %>%
        group_by(CDR3.nt, CDR3.aa, V.name, D.name, J.name, Sequence, chain) %>%
        summarise(
            Clones=sum(Clones),
            Proportion=sum(Proportion),
            V.end=V.end[1],
            D.start=D.start[1],
            D.end =D.end[1],
            J.start=J.start[1],
            VJ.ins=VJ.ins[1],
            VD.ins=VD.ins[1],
            DJ.ins=DJ.ins[1],
            Barcode=paste(Barcode, collapse=";"),
            raw_clonotype_id=paste(raw_clonotype_id, collapse=";"),
            ContigID=paste(ContigID, collapse=";"),
        ) %>%
        as.data.frame()
}

immdata$meta = left_join(
    immdata$meta,
    metadata %>% select(-"TCRDir"),
    by = "Sample"
)

saveRDS(immdata, file=rdsfile)
