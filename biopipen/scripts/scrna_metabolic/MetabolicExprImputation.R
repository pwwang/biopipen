library(scImpute)
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)

srtdir = {{in.srtdir | r}}
outdir = "{{out.outdir}}/"
joboutdir = {{job.outdir | r}}
infmt = {{envs.infmt | r}}
drop_thre = {{envs.drop_thre | r}}

ncores = {{envs.ncores | r}}

do_one_srtfile = function(srtfile) {
    prefix = tools::file_path_sans_ext(basename(srtfile))
    outfile = file.path(outdir, paste0(prefix, ".RDS"))
    cachedfile = file.path(joboutdir, paste0(prefix, ".RDS"))

    if (file.exists(cachedfile)) {
        file.symlink(normalizePath(cachedfile), outfile)
        return(NULL)
    }
    sobj = readRDS(srtfile)
    counts = as.data.frame(sobj@assays$RNA@counts)
    labels = as.integer(Idents(sobj))

    count_path = file.path(
        outdir,
        paste0(prefix, "_rawcounts.rds")
    )
    saveRDS(counts, file=count_path)

    if (any(table(labels)) < 10) {
        scimpute(
            count_path,
            infile = "rds",
            outfile = "rds",
            type = "count",
            out_dir = outdir,
            labeled = FALSE,
            drop_thre = drop_thre,
            labels = NULL,
            Kcluster = length(unique(labels)),
            ncores = ncores
        )
    } else {
        scimpute(
            count_path,
            infile = "rds",
            outfile = "rds",
            type = "count",
            out_dir = outdir,
            labeled = TRUE,
            drop_thre = drop_thre,
            labels = labels,
            ncores = ncores
        )
    }

    ofile = file.path(outdir, "scimpute_count.rds")
    file.rename(ofile, cachedfile)
    file.symlink(normalizePath(cachedfile), outfile)
}

for (srtfile in Sys.glob(file.path(srtdir, "*.RDS"))) {
    do_one_srtfile(srtfile)
}
