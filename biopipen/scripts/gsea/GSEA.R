# devtools::install_github("GSEA-MSigDB/GSEA_R")

{{ biopipen_dir | joinpaths: "utils", "io.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "gsea.R" | source_r }}

library(dplyr)
library(tibble)
library(GSEA)

# ?GSEA to see help

infile <- {{ in.infile | r }}
metafile <- {{ in.metafile | r }}
gmtfile <- {{ in.gmtfile | r }}
{% if in.configfile %}
config = {{in.config | read | toml_loads | r}}
{% else %}
config = list()
{% endif %}
outdir <- {{ out.outdir | r }}
envs <- {{ envs | r }}

clscol <- if (is.null(config$clscol)) envs$clscol else config$clscol
inopts <- envs$inopts
envs$inopts <- NULL
metaopts <- envs$metaopts
envs$metaopts <- NULL
envs$clscol <- NULL
if (!is.null(config$doc_string)) {
    envs$doc.string = config$doc_string
}
envs$doc_string <- NULL

if (is.character(inopts) && inopts == "rds") {
    indata = readRDS(infile)
} else {
    indata = read.table.opts(infile, inopts)
}

metadata = read.table.opts(metafile, metaopts)
classes = metadata[colnames(indata), clscol]

runGSEA(
    indata, # expression data
    classes, # sample classes
    gmtfile, # the GMT file
    outdir,
    envs
)
