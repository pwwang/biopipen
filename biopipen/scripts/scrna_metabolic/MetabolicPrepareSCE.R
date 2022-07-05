source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/rnaseq.R")

library(scater)
library(Seurat)

impfiles = {{in.impfiles | r}}
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

sceobjs = list()
by.rmagic = FALSE
for (i in seq_len(length(impfiles))) {
    sobj = readRDS(impfiles[i])
    if ("UNIMPUTED_RNA" %in% Assays(sobj)) {
        by.rmagic = TRUE
    }
    if (by.rmagic) {
        tpms = GetAssayData(sobj, assay = "RNA", slot = "data")
    } else {
        counts = GetAssayData(sobj, "counts")
        tpms = unit_conversion(
            counts,
            "count",
            "tpm",
            args = list(genelen = glenFromGFFExons(refexon))
        )
        # counts = counts[rownames(tpms),]
    }
    srt =  as.SingleCellExperiment(sobj, assay="RNA")
    cdata = colData(srt)
    bname = tools::file_path_sans_ext(basename(impfiles[i]))
    # jobname.case_subset
    cdata$.subset = sub("^.+?_", "", bname)
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

sce = do_call(cbind, sceobjs)
if (by.rmagic) {
    attr(sce, "by.rmagic") = TRUE
}
saveRDS(sce, file=outfile)
