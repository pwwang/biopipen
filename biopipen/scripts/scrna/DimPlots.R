library(Seurat)

srtfile = {{in.srtobj | quote}}
groupfile = {{in.groupfile | r}}
config = {{in.configfile | read | toml_loads | r}}
outfile = {{out.outfile | quote}}

set.seed(8525)

sobj = readRDS(srtfile)
config$name = NULL

if (config$by.ident) {
    # by ident
    config$by.ident = NULL
    config$object = sobj
    config$group.by = "ident"
    p = do.call(DimPlot, config)
    png(outfile, res=100, width=1000, height=1000)
    print(p)
    dev.off()
} else {
    # TODO
}
