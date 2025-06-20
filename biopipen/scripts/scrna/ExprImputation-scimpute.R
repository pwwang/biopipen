library(scImpute)
library(Seurat)

infile = {{in.infile | r}}
outfile = {{out.outfile | r}}
joboutdir = {{job.outdir | append: "/" | r}}
drop_thre = {{envs.scimpute_args.drop_thre | r}}
kcluster = {{(envs.scimpute_args.kcluster | default: None | r}}
ncores = {{envs.scimpute_args.ncores | r}}
refgene = {{envs.scimpute_args.refgene | r}}

setwd(joboutdir)

labels = NULL
sobj = read_obj(infile)
counts = as.data.frame(sobj@assays$RNA@counts)
kc = length(unique(Idents(sobj)))
if (kc > 0) {
    labels = as.integer(Idents(sobj))
}

count_path = file.path(joboutdir, "counts.rds")
saveRDS(counts, file=count_path)

scimpute(
    count_path,
    infile = "rds",
    outfile = "rds",
    type = "count",
    out_dir = joboutdir,
    labeled = !is.null(labels),
    drop_thre = drop_thre,
    Kcluster = kcluster,
    labels = labels,
    ncores = ncores,
)

imputed = readRDS(file.path(joboutdir, "scimpute_count.rds"))
outobj = CreateSeuratObject(counts = imputed)

outobj@meta.data = sobj@meta.data[rownames(outobj@meta.data), , drop=FALSE]
# remember that it is the counts being imputed, we still need to
# normalize the data
outobj@misc$impute_method = "scimpute"

save_obj(outobj, outfile)
