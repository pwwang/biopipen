{{rimport}}('__init__.r')

infiles  = {{i.infiles | R}}
outfile  = {{o.outfile | R}}
inopts   = {{args.inopts, i.infiles | lambda x, files: (x, len(files)) | lambda x: {k:(v if isinstance(v, list) and len(v) > 1 else v*x[1] if isinstance(x, list) and len(v) == 1 else [v]*x[1]) for k, v in x[0].items()} | R}}
params   = {{args.params | R}}
na       = {{args.na | R}}
fn2rname = {{args.fn2rname}}
fill     = as.logical({{args.fill | R}})

mats = list()
for (i in 1:length(infiles)) {
	inparams = c(list(
		file      = infiles[i],
		header    = ifelse(is.null(inopts$cnames[i]), T, as.logical(inopts$cnames[i])),
		row.names = ifelse(is.null(inopts$rnames[i]) || inopts$rnames[i], 1, NULL),
		sep       = ifelse(is.null(inopts$delimit[i]), "\t", inopts$delimit[i]),
		skip      = ifelse(is.null(inopts$skip[i]), 0, inopts$skip[i])
	), params)
	mats[[i]]  = do.call(read.table, inparams)

	if (is.null(inparams$row.names) && !is.null(fn2rname)) {
		rname = fn2rname(tools::file_path_sans_ext(basename(infiles[i])))
		if (nrow(mats[[i]]) == 1)
			rownames(mats[[i]]) = rname
		else
			rownames(mats[[i]]) = paste(rname, 1:nrow(mats[[i]]), sep='_r')
	}
}
mat = ifelse(fill, do.call(rbind.fill, mats), do.call(rbind, mats))
mat[is.na(mat)] = na

outparams = list(
	x         = mat,
	file      = outfile,
	sep       = '\t',
	quote     = F,
	col.names = ifelse(is.null(inopts$cnames[1]), T, as.logical(inopts$cnames[1])),
	row.names = T
)
do.call(write.table, outparams)
