# devtools::install_github("GSEA-MSigDB/GSEA_R")

source("{{biopipen_dir}}/utils/io.R")
source("{{biopipen_dir}}/utils/gsea.R")

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
if (!is.null(config$doc.string)) {
    envs$doc.string = config$doc.string
}

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
