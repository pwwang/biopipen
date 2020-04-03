{{"__init__.R" | rimport}}

infile  = {{i.infile | R}}
outdir  = {{o.outdir | R}}
inopts  = {{args.inopts | R}}
size    = {{args.size | int | R}}

mat = read.table.inopts(infile, inopts)

allrnames = if (as.logical(inopts$rnames)) rownames(mat) else paste0("ROW", 1:nrow(mat))

for ( chunk in split(1:length(allrnames), ceiling(seq_along(1:length(allrnames))/size)) ) {
	rnames    = allrnames[chunk]
	m         = mat[rnames, , drop = F]
	fn        = paste(gsub("[[:punct:]]", "_", rnames), collapse='-')
	fn        = paste0({{i.infile | fn2 | quote}}, '_', fn, {{i.infile | ext | quote}})
	outfile   = file.path(outdir, fn)
	outparams = list(
		x         = m,
		file      = outfile,
		sep       = list.get(inopts, "delimit", "\t"),
		quote     = F,
		col.names = is.true(inopts$cnames),
		row.names = is.true(inopts$rnames)
	)
	do.call(write.table, outparams)
}
