source("{{biopipen_dir}}/utils/misc.R")

library(SeuratWrappers)
library(Seurat)

infile = {{in.infile | r}}
outfile = {{out.outfile | r}}
envs = {{envs.alra_args | r}}

log_info("Loading Seurat object")
sobj = readRDS(infile)
assay <- DefaultAssay(sobj)

log_info("Imputing expression values, using ALRA")
envs$object = sobj
sobj = do_call(RunALRA, envs)

# sobj = RunALRA(sobj)
log_info("Renaming assays")
sobj = RenameAssays(sobj, assay.name = assay, new.assay.name = "RAW")
sobj = RenameAssays(sobj, assay.name = "alra", new.assay.name = assay)
DefaultAssay(sobj) <- assay

sobj@misc$impute = "alra"

log_info("Saving Seurat object")
saveRDS(sobj, outfile)

# choosek_plot_file = file.path(dirname(outfile), "choosek.png")
# png(choosek_plot_file, width = 1200, height = 1000, res = 100)
# p = ALRAChooseKPlot(sobj)
# print(p)
# dev.off()
