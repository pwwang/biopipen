library(scImpute)
library(dplyr)
library(tidyr)
library(tibble)
library(rtracklayer)

infile = {{in.infile | r}}
groupfile = {{in.groupfile | r}}
outfile = {{out.outfile | r}}
joboutdir = "{{job.outdir}}/"
infmt = {{envs.infmt | r}}
outfmt = {{envs.outfmt | r}}
datatype = {{envs.type | r}}
drop_thre = {{envs.drop_thre | r}}
kcluster = {{(envs.kcluster or envs.Kcluster or None) | r}}

ncores = {{envs.ncores | r}}
refgene = {{envs.refgene | r}}

labels = NULL
if (infmt == "seurat") {
    library(Seurat)
    sobj = readRDS(infile)
    counts = as.data.frame(sobj@assays$RNA@counts)
    datatype = "count"
    kc = length(unique(Idents(sobj)))
    if (kc > 0) {
        labels = as.integer(Idents(sobj))
    }
} else {
    if (infmt == "rds") {
        counts = readRDS(infile)
    } else if (infmt == "txt") {
        counts = read.table(
            infile,
            header=T, sep="\t", row.names=1, quote=F, check.names = F
        )
    } else if (infmt == "csv") {
        counts = read.table(
            infile,
            header=T, sep=",", row.names=1, quote=F, check.names = F
        )
    } else {
        stop(paste("Unsupported `envs.infmt`", infmt))
    }
}

if (!is.null(groupfile)) {
    groups = read.table(groupfile, row.names=NULL, header=T, sep="\t", check.names = F)
    n_samples = ncol(groups) - 1
    groupcol = colnames(groups)[1]
    labels = groups %>% rowwise() %>%
        mutate(
            across(2:(n_samples+1),
            ~ sapply(
                strsplit(.x, ";", fixed=TRUE),
                function(x) paste(cur_column(), x, sep="_", collapse=";")
            )
        )) %>%
        unite("ALL", 2:(n_samples+1), sep=";") %>%
        # group_by(across(all_of(groupcol))) %>%
        # summarise(across(everything(), ~ paste(.x, collapse=";")))
        rowwise() %>%
        group_map(
            function(x, ...) {
                cells = unlist(strsplit(x$ALL, ";"))
                out = rep(x$Group, length(cells))
                names(out) = cells
                out
            }
        ) %>% unlist()

    cells = intersect(names(labels), colnames(counts))
    labels = labels[cells]
    labels = as.integer(as.factor(labels))
    counts = counts[, cells, drop=F]
}

if (datatype != "count") {
    genelen_df = as.data.frame(rtracklayer::import(refgene)) %>%
        column_to_rownames("gene_id") %>%
        select("width")
    genes = intersect(rownames(counts), rownames(genelen_df))
    genelen = genelen_df[genes, "width"]
} else {
    genelen = NULL
}

count_path = file.path(joboutdir, "counts.rds")
saveRDS(counts, file=count_path)

scimpute(
    count_path,
    infile = "rds",
    outfile = outfmt,
    type = datatype,
    out_dir = joboutdir,
    labeled = !is.null(labels),
    drop_thre = drop_thre,
    Kcluster = kcluster,
    labels = labels,
    ncores = ncores,
    genelen = genelen
)

ofile = file.path(joboutdir, paste0("scimpute_count.", outfmt))
file.rename(ofile, outfile)
