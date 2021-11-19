# devtools::install_github("GSEA-MSigDB/GSEA_R")
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

replace_underscore_keys = function(alist) {
    newlist = alist
    names(newlist) = sub("_", ".", names(alist), fixed=T)
    return(newlist)
}

config = replace_underscore_keys(config)
envs = replace_underscore_keys(envs)
clscol <- if (is.null(config$clscol)) envs$clscol else config$clscol
infmt <- envs$infmt
envs$infmt <- NULL
envs$clscol <- NULL
if (!is.null(config$doc.string)) {
    envs$doc.string = config$doc.string
}

#
runGSEA = function(
    # The expression file
    infile,
    # input format: matrix or rds
    infmt,
    # The metadata file
    metafile,
    # the column from metafile for classes
    clscol,
    # The gmt file
    gmtfile,
    # arguments for GSEA()
    envs,
    # output directory
    outdir
) {

    if (is.null(clscol)) {
        stop("`envs.clscol` is required to determine CLS from in.metafile")
    }

    # reproducibility
    if (is.null(envs$random.seed)) {
        envs$random.seed <- 8525
    }

    # prepare gct file
    gctfile = file.path(outdir, "gsea.gct")
    if (infmt == "matrix") {
        if (endsWith(infile, ".gz")) {
            indata <- read.table(
                gzfile(infile),
                header = T, row.names = NULL, sep = "\t", check.names = F
            )
        } else {
            indata <- read.table(
                infile,
                header = T, row.names = NULL, sep = "\t", check.names = F
            )
        }
        genes = make.unique(indata[, 1])
        indata = indata[, -1, drop=F]
        rownames(indata) = genes
    } else {
        indata <- readRDS(infile)
    }
    con = file(gctfile, open='w')
    write("#1.2", con)
    write(paste(dim(indata), collapse = "\t"), con)
    close(con)
    indata = indata %>%
        mutate(Description = "na") %>%
        rownames_to_column("NAME") %>%
        select(NAME, Description, everything())
    write.table(
        indata,
        gctfile,
        row.names = F,
        col.names = T,
        sep="\t",
        quote=F,
        append = T
    )

    # prepare cls file
    clsfile = file.path(outdir, "gsea.cls")
    metadata <- read.table(
        metafile,
        header = T,
        row.names = NULL,
        sep = "\t",
        check.names = F
    ) %>% select(Sample, all_of(clscol)) %>% column_to_rownames("Sample")
    allclasses = as.character(
        metadata[colnames(indata)[3:ncol(indata)], 1, drop=T]
    )
    classes = unique(allclasses)
    con = file(clsfile, open='w')
    write(paste(length(allclasses), length(classes), '1'), con)
    write(paste('#', paste(classes, collapse=" ")), con)
    write(paste(allclasses, collapse=" "), con)
    close(con)

    envs$input.ds = gctfile
    envs$input.cls = clsfile
    envs$gs.db = gmtfile
    envs$output.directory = outdir

    do.call(GSEA, envs)
}

runGSEA(
    infile=infile,
    infmt=infmt,
    metafile=metafile,
    clscol=clscol,
    gmtfile=gmtfile,
    envs=envs,
    outdir=outdir
)
