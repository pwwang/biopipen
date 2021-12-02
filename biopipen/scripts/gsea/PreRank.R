# PreRank the genes for GSEA analysis
# See: https://gseapy.readthedocs.io/en/latest/_modules/gseapy/algorithm.html#ranking_metric
source("{{biopipen_dir}}/utils/io.R")
source("{{biopipen_dir}}/utils/gsea.R")

infile = {{in.infile | quote}}
metafile = {{in.metafile | quote}}
{% if in.configfile %}
config = {{in.config | read | toml_loads | r}}
{% else %}
config = list()
{% endif %}
outfile = {{out.outfile | quote}}
envs = {{envs | r}}
clscol <- if (is.null(config$clscol)) envs$clscol else config$clscol
classes <- if (is.null(config$classes)) envs$classes else config$classes

if (is.null(clscol)) {
    stop("No `clscol` specified.")
}

if (is.null(classes) || length(classes) %% 2 != 0) {
    stop(paste("`classes` must be pair(s) of labels."))
}

if (is.character(envs$inopts) && inopts == "rds") {
    indata = readRDS(infile)
} else {
    indata = read.table.opts(infile, envs$inopts)
}

metadata = read.table.opts(metafile, envs$metaopts)
allclasses = metadata[colnames(indata), clscol]

out = NULL
for (labels in split(classes, ceiling(seq_along(classes)/2))) {
    pos = labels[1]
    neg = labels[2]
    rnk = prerank(indata, pos, neg, allclasses, envs$method)
    if (is.null(out)) {
        out = rnk
    } else {
        out = full_join(out, rnk, by="Gene")
    }
}

write.table(out, outfile, row.names=F, col.names=T, sep="\t", quote=F)

