{{ biopipen_dir | joinpaths: "utils", "io.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "plot.R" | source_r }}

infile = {{in.infile | r}}
outfile = {{out.outfile | r}}
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
    indata = apply(indata, 2, function(x) elems[as.logical(x)])
}

plotVenn(
    indata,
    args,
    ggs,
    devpars,
    outfile
)
