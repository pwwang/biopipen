source("{{biopipen_dir}}/utils/misc.R")

library(SeuratWrappers)
library(Seurat)

infile = {{in.infile | r}}
outfile = {{out.outfile | r}}
envs = {{envs.alra_args | r}}

print("Loading Seurat object")
sobj = readRDS(infile)
DefaultAssay(sobj) <- "RNA"

print("Imputing expression values, using ALRA")
envs$object = sobj
sobj = do_call(RunALRA, envs)

# sobj = RunALRA(sobj)
print("Renaming assays")
sobj = RenameAssays(sobj, RNA = "UNIMPUTED_RNA")
sobj = RenameAssays(sobj, alra = "RNA")
DefaultAssay(sobj) <- "RNA"

attr(sobj, "impute") = "alra"
print("Saving Seurat object")
saveRDS(sobj, outfile)

# choosek_plot_file = file.path(dirname(outfile), "choosek.png")
# png(choosek_plot_file, width = 1200, height = 1000, res = 100)
# p = ALRAChooseKPlot(sobj)
# print(p)
# dev.off()
