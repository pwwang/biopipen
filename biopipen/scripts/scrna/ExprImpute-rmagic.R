
tryCatch({
    # in order to load Rmagic
    workdir = {{job.outdir | r}}
    conda_prefix = Sys.getenv("CONDA_PREFIX")
    setwd(workdir)
    file.symlink(conda_prefix, "miniconda3")
}, error=function(e) {})

Sys.setenv(RETICULATE_PYTHON = {{envs.rmagic_args.python | r}})
reticulate::use_python({{envs.rmagic_args.python | r}})

library(Rmagic)
library(Seurat)

infile = {{in.infile | r}}
outfile = {{out.outfile | r}}

sobj = readRDS(infile)
DefaultAssay(sobj) <- "RNA"

sobj = magic(sobj)
sobj = RenameAssays(sobj, RNA = "RAW_RNA", MAGIC_RNA = "RNA")

saveRDS(sobj, outfile)
