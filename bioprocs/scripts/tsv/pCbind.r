{{rimport}}('__init__.r')

infiles  = {{in.infiles | R}}
outfile  = {{out.outfile | R}}
inopts   = {{args.inopts, in.infiles | lambda x, files: (x, len(files)) | lambda x: {k:(v if isinstance(v, list) and len(v) > 1 else v*x[1] if isinstance(x, list) and len(v) == 1 else [v]*x[1]) for k, v in x[0].items()} | R}}
params   = {{args.params | R}}
na       = {{args.na | R}}
fn2cname = {{args.fn2cname}}
fill     = as.logical({{args.fill | R}})

bindfunc = if(fill) cbind.fill else cbind

mat = NULL
for (i in 1:length(infiles)) {
	inparams = c(list(
		file      = infiles[i],
		header    = if(is.null(inopts$cnames[i])) T else as.logical(inopts$cnames[i]),
		row.names = if(is.null(inopts$cnames[i]) || inopts$cnames[i]) 1 else NULL,
		sep       = if(is.null(inopts$delimit[i])) "\t" else inopts$delimit[i]
	), params)
	.mat  = do.call(read.table, inparams)
	if (!inparams$header)
		colnames(.mat) = fn2cname(tools::file_path_sans_ext(infiles[i]))
	mat   = if(is.null(mat)) .mat else bindfunc(mat, .mat)
}

mat[is.na(mat)] = na

outparams = list(
	x         = mat,
	file      = outfile,
	sep       = if(is.null(inopts$delimit[1])) "\t" else inopts$delimit[1],
	quote     = F,
	col.names = T,
	row.names = if(is.null(inopts$rnames[1])) T else as.logical(inopts$rnames[1])
)
do.call(write.table, outparams)
