source("{{biopipen_dir}}/utils/rnaseq.R")

library(scater)
library(Seurat)

impdir = {{in.impdir | r}}
srtdir = {{in.srtdir | r}}
gmtfile = {{in.gmtfile | r}}
refexon = {{envs.refexon | r}}
outfile = {{out.outfile | r}}

## gmtPathways is copied from fgsea package.
gmtPathways <- function(gmt.file) {
    pathwayLines <- strsplit(readLines(gmt.file), "\t")
    pathways <- lapply(pathwayLines, tail, -2)
    names(pathways) <- sapply(pathwayLines, head, 1)
    pathways
}
pathways = gmtPathways(gmtfile)
metabolics = unique(as.vector(unname(unlist(pathways))))

impfiles = sort(Sys.glob(file.path(impdir, "*.RDS")))
srtobjfiles = sort(Sys.glob(file.path(srtdir, "*.RDS")))

sceobjs = list()
for (i in seq_len(length(impfiles))) {
    counts = readRDS(impfiles[i])
    tpms = unit_conversion(
        counts,
        "count",
        "tpm",
        args = list(genelen = glenFromGFFExons(refexon))
    )
    counts = counts[rownames(tpms),]
    srt =  as.SingleCellExperiment(readRDS(srtobjfiles[i]), assay="RNA")
    cdata = colData(srt)
    cdata$.subset = tools::file_path_sans_ext(basename(impfiles[i]))
    rdata = rowData(srt)[rownames(tpms),,drop=F]
    rdata$metabolic = FALSE
    rdata[rownames(rdata) %in% metabolics, "metabolic"] = TRUE
    sceobj = SingleCellExperiment(
        assays = list(tpm=data.matrix(tpms), exprs=data.matrix(log2(tpms + 1))),
        colData = cdata,
        rowData = rdata
    )
    sceobjs[[i]] = sceobj
}

sce = do.call(cbind, sceobjs)
saveRDS(sce, file=outfile)
