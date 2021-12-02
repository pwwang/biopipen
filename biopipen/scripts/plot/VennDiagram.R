library(dplyr)

source("{{biopipen_dir}}/utils/io.R")
source("{{biopipen_dir}}/utils/plot.R")

infile = {{in.infile | quote}}
outfile = {{out.outfile | quote}}
inopts = {{envs.inopts | r}}
intype = {{envs.intype | r}}
devpars = {{envs.devpars | r}}
args = {{envs.args | r}}
ggs = {{envs.ggs | r}}

indata = read.table.opts(infile, inopts)

if (intype == "raw") {
    indata = setNames(indata[[1]], rownames(indata))
    indata = lapply(indata, function(x) unlist(strsplit(x, ",", fixed=TRUE)))
} else { # computed
    elems = rownames(indata)
    indata = indata %>%
        mutate(across(everything(), function(x) elems[as.logical(x)])) %>%
        as.list()
}

plotVenn(
    indata,
    args,
    ggs,
    devpars,
    outfile
)
