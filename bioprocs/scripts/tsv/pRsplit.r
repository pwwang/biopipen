{{rimport}}('__init__.r')

infile  = {{i.infile | R}}
outdir  = {{o.outdir | R}}
params  = {{args.params | R}}
inopts  = {{args.inopts | R}}
size    = as.numeric({{args.size | R}})

inparams = list(
	file      = infile,
	header    = as.logical(inopts$cnames),
	row.names = if (as.logical(inopts$rnames)) 1 else NULL,
	skip      = if (is.null(inopts$skip)) 0 else as.numeric(inopts$skip),
	sep       = inopts$delimit
)

mat       = do.call(read.table.nodup, c(inparams, params))
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
		sep       = inopts$delimit,
		quote     = F,
		col.names = as.logical(inopts$cnames),
		row.names = as.logical(inopts$rnames)
	)
	do.call(write.table, outparams)
}
