library(SeuratWrappers)
library(Seurat)

infile = {{in.infile | r}}
outfile = {{out.outfile | r}}
envs = {{envs.alra_args | r}}

sobj = readRDS(infile)
DefaultAssay(sobj) <- "RNA"

envs$object = sobj
sobj = do.call(RunALRA, envs)

sobj = RunALRA(sobj)
sobj = RenameAssays(sobj, RNA = "UNIMPUTED_RNA")
sobj = RenameAssays(sobj, alra = "RNA")
DefaultAssay(sobj) <- "RNA"

attr(sobj, "impute") = "alra"
saveRDS(sobj, outfile)

# choosek_plot_file = file.path(dirname(outfile), "choosek.png")
# png(choosek_plot_file, width = 1200, height = 1000, res = 100)
# p = ALRAChooseKPlot(sobj)
# print(p)
# dev.off()
