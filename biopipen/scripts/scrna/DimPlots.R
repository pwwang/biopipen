library(Seurat)
library(dplyr)

srtfile = {{in.srtobj | r}}
{% if in.configfile %}
config = {{in.configfile | read | toml_loads | r}}
{% else %}
config = list()
{% endif %}
outdir = {{out.outdir | r}}
cases = {{envs.cases | r}}

set.seed(8525)

sobj = readRDS(srtfile)

for (case in names(config$cases)) {
    cases[[case]] = config$cases[[case]]
}

for (case in names(cases)) {
    args = cases[[case]]
    if (!is.null(args$mutate)) {
        if (args$group.by %in% colnames(sobj@meta.data)) {
            i = 1
            group.by = paste0(args$group.by, "..", i)
            while (group.by %in% colnames(sobj@meta.data)) {
                i = i + 1
                group.by = paste0(args$group.by, "..", i)
            }
            args$group.by = group.by
        }
        sobj@meta.data = sobj@meta.data %>% mutate(
            ..new.meta = eval(parse(text=args$mutate))
        )
        colnames(sobj@meta.data)[ncol(sobj@meta.data)] = args$group.by
        args$mutate = NULL
    }
    if (!is.null(args$reduction) && !args$reduction %in% names(sobj@reductions)) {
        if (args$reduction == "umap") {
            sobj = RunUMAP(sobj)
        } else if (args$reduction == "tsne") {
            sobj = RunTSNE(sobj)
        } else if (args$reduction == "pca") {
            sobj = RunPCA(sobj)
        }
    }
    args$object = sobj
    p = do.call(DimPlot, args)
    outfile = file.path(outdir, paste0(case, ".png"))
    png(outfile, res=100, width=1000, height=1000)
    print(p)
    dev.off()
}
