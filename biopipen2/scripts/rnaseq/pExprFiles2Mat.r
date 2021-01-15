{{rimport}}('__init__.r')
options(stringsAsFactors = FALSE)

infiles   = {{ i.infiles | R}}
outfile   = {{ o.outfile | quote}}
inopts    = {{ args.inopts | R}}
fn2sample = {{args.fn2sample}}

ret = NULL
for (infile in infiles) {
	sample = fn2sample(tools::file_path_sans_ext(basename(infile)))
	indata = read.table.inopts(infile, inopts)[, 1, drop = FALSE]
	colnames(indata) = sample
	if (is.null(ret)) {
		ret = indata
	} else {
		ret = cbind(ret, indata[rownames(ret),,drop = FALSE])
	}
}
write.table(ret, outfile, sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
