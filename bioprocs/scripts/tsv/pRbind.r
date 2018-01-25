{{rbindfill}}

cnames = as.logical({{args.cnames | lambda x: x if isinstance(x, list) else [x] | Rvec}})
rnames = as.logical({{args.rnames | lambda x: x if isinstance(x, list) else [x] | Rvec}})

mat = NULL
infiles = {{in.infiles | Rvec}}
inlen   = length(infiles)
cnames  = if (length(cnames) == 1) rep(cnames, inlen) else c(cnames, rep(TRUE, inlen - length(cnames)))
rnames  = if (length(rnames) == 1) rep(rnames, inlen) else c(rnames, rep(TRUE, inlen - length(rnames)))
for (i in 1:length(infiles)) {
	infile = infiles[i]
	mat2 = read.table (infile, sep="\t", header = cnames[i], row.names = if (rnames[i]) 1 else NULL, check.names = F)
	mat  = rbindfill(mat, mat2)
}
mat[is.na(mat)] = {{args.na | R}}

write.table (mat, 
	{{out.outfile | quote}}, 
	sep="\t", 
	quote=F, 
	col.names = as.logical({{args.cnames | lambda x: x if isinstance(x, list) else [x] | lambda x: any(x) | R}}), 
	row.names = as.logical({{args.rnames | lambda x: x if isinstance(x, list) else [x] | lambda x: any(x) | R}})
)