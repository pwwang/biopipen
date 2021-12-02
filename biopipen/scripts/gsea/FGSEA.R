# PreRank the genes for GSEA analysis
# See: https://gseapy.readthedocs.io/en/latest/_modules/gseapy/algorithm.html#ranking_metric
source("{{biopipen_dir}}/utils/io.R")
source("{{biopipen_dir}}/utils/gsea.R")

infile = {{in.infile | quote}}
metafile = {{in.metafile | quote}}
gmtfile = {{in.gmtfile | quote}}
{% if in.configfile %}
config = {{in.config | read | toml_loads | r}}
{% else %}
config = list()
{% endif %}
outdir = {{out.outdir | quote}}
envs = {{envs | r}}
clscol <- if (is.null(config$clscol)) envs$clscol else config$clscol
classes <- if (is.null(config$classes)) envs$classes else config$classes

if (is.null(clscol)) {
    stop("No `clscol` specified.")
}

if (is.null(classes) || length(classes) != 2) {
    stop(paste("`classes` must be a pair of labels."))
}

if (is.character(envs$inopts) && inopts == "rds") {
    indata = readRDS(infile)
} else {
    indata = read.table.opts(infile, envs$inopts)
}

metadata = read.table.opts(metafile, envs$metaopts)
allclasses = metadata[colnames(indata), clscol]

ranks = prerank(indata, classes[1], classes[2], allclasses, envs$method)

write.table(
    ranks,
    file.path(outdir, "fgsea.rank"),
    row.names=F,
    col.names=T,
    sep="\t",
    quote=F
)

top = envs$top
envs$nproc = envs$ncores
envs$inopts = NULL
envs$metaopts = NULL
envs$method = NULL
envs$clscol = NULL
envs$classes = NULL
envs$ncores = NULL
envs$top = NULL
# the rest are the arguments for `fgsea()`

runFGSEA(ranks, gmtfile, top, outdir, envs)
