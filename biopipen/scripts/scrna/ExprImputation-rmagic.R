
tryCatch({
    # in order to load Rmagic
    workdir = {{job.outdir | r}}
    conda_prefix = Sys.getenv("CONDA_PREFIX")
    setwd(workdir)
    file.symlink(conda_prefix, "miniconda3")
}, error=function(e) {})

python = {{envs.rmagic_args.python | r}}
Sys.setenv(RETICULATE_PYTHON = Sys.which(python))
# reticulate::use_python(python)

library(Rmagic)
library(Seurat)

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}

sobj <- readRDS(infile)
assay <- DefaultAssay(sobj)

sobj <- magic(sobj)
sobj <- RenameAssays(sobj, assay.name = assay, new.assay.name = "RAW")
sobj <- RenameAssays(sobj, assay.name = "MAGIC_RNA", new.assay.name = assay)

DefaultAssay(sobj) <- assay

attr(sobj, "impute") = "rmagic"
saveRDS(sobj, outfile)
