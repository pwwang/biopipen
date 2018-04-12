{{rimport}}('__init__.r')


inopts.cnames = {{args.inopts.cnames | lambda x: x if isinstance(x, list) else [x] | Rvec}}
inopts.rnames = {{args.inopts.rnames | lambda x: x if isinstance(x, list) else [x] | Rvec}}
inopts.cnames = as.logical(inopts.cnames)
inopts.rnames = as.logical(inopts.rnames)

infiles = {{in.infiles | Rvec}}
outfile = {{out.outfile | quote}}
na      = {{args.na | R}}

mat = NULL
inlen   = length(infiles)
cnames  = if (length(inopts.cnames) == 1) rep(inopts.cnames, inlen) else c(inopts.cnames, rep(TRUE, inlen - length(inopts.cnames)))
rnames  = if (length(inopts.rnames) == 1) rep(inopts.rnames, inlen) else c(inopts.rnames, rep(TRUE, inlen - length(inopts.rnames)))
for (i in 1:inlen) {
	infile = infiles[i]
	mat2 = read.table (infile, sep="\t", header = cnames[i], row.names = if (rnames[i]) 1 else NULL, check.names = F)
	mat  = rbindfill(mat, mat2)
}
mat[is.na(mat)] = na

write.table (mat, 
	outfile, 
	sep="\t", 
	quote=F, 
	col.names = any(cnames), 
	row.names = any(rnames)
)