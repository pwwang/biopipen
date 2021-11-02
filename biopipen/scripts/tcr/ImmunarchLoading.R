# Loading 10x data into immunarch
library(immunarch)
library(dplyr)
library(tidyr)

metafile = {{ in.metafile | quote }}
rdsfile = {{ out.rdsfile | quote }}
allclfile = {{ out.allclones | quote }}
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

datadir = tempfile(pattern = "immunarch-", tmpdir = tmpdir)
dir.create(datadir, showWarnings = FALSE)

# Find filtered_contig_annotations.csv and link then in datadir
for (i in seq_len(nrow(metadata))) {
    sample = as.character(metadata[i, "Sample"])
    annofile = file.path(
        as.character(metadata[i, "TCRDir"]),
        "filtered_contig_annotations.csv"
    )
    if (!file.exists(annofile)) {
        stop(paste(
            "Cannot find `filtered_contig_annotations.csv`",
            "in given TCRDir for sample:",
            sample
        ))
    }

    file.symlink(
        normalizePath(annofile),
        file.path(datadir, paste0(sample, ".csv"))
    )
}

immdata = repLoad(datadir)
# drop TRAs for paired data
immdata$single = list()
for (sample in names(immdata$data)) {
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

cdr3aa = c()
for (sample in names(immdata$data)) {
    cdr3aa = unique(c(cdr3aa, immdata$data[[sample]]$CDR3.aa))
}
out = data.frame(CDR3.aa = cdr3aa)
for (sample in names(immdata$data)) {
    out = left_join(out, immdata$data[[sample]], by="CDR3.aa") %>%
        select(all_of(colnames(out)), Barcode)
    colnames(out)[ncol(out)] = sample
}
out = out %>% distinct(CDR3.aa, .keep_all = TRUE)

write.table(out, allclfile, col.names=T, row.names=F, sep="\t", quote=F)
