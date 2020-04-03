{{"__init__.R", "plot.R" | rimport}}
options(stringsAsFactors = FALSE)

infile  = {{ i.infile | quote}}
outfile = {{ o.outfile | quote}}
hifile1 = {{ i.hifile | R}}
hifile2 = {{ args.hifile | R}}
inopts  = {{ args.inopts | R}}
devpars = {{ args.devpars | R}}
# allow bin_size instead of bin.size
params  = {{ args.params | .items:
                         | :[(item[0].replace("_", "."), item[1])
                             for item in _]
                         | dict | R}}

indata = read.table.inopts(infile, inopts)

hifile = NULL
if (is.true(hifile1)) {
    hifile = hifile1
} else if (is.true(hifile2)) {
    hifile = hifile2
}
if (is.true(hifile) && is.false(params$highlight)) {
    hidata = read.table.inopts(hifile, list(rnames=FALSE, cnames=FALSE))
    hidata = hidata[ hidata[,1] %in% indata[,1], , drop=FALSE]
    params$highlight = hidata[,1]
    if (ncol(hidata) > 1) {
        params$highlight.text = hidata[,2]
    }
}


plot.man2(
	indata,
	plotfile = outfile,
	params   = params,
	devpars  = devpars
)
