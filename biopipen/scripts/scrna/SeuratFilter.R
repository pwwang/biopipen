library(Seurat)
library(rlang)
library(dplyr)

srtobjfile = {{in.srtobj | r}}
outfile = {{out.outfile | r}}
mutaters = {{in.filters | config: "toml" | attr: "mutaters" | r}}

srtobj = readRDS(srtobjfile)


if (!is.null(mutaters)) {
    expr = list()
    for (key in names(mutaters)) {
        expr[[key]] = parse_expr(mutaters[[key]])
    }
    srtobj@meta.data = srtobj@meta.data |> mutate(!!!expr)
}


sobj = subset(srtobj, subset = {{in.filters | config: "toml" | attr: "filter"}})
saveRDS(sobj, outfile)
