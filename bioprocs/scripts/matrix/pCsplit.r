cnames = as.logical({{args.cnames | R}})
rnames = as.logical({{args.rnames | R}})

mat    = read.table ({{in.infile | quote}}, sep="\t", header = cnames, row.names = if (rnames) 1 else NULL, check.names = F)
heads  = if (cnames) colnames(mat) else paste("COL", 1:ncol(mat), sep="")
for (i in 1:ncol(mat)) {
	cname = heads[i]
	m     = mat[, cname, drop=F]
	ofile = file.path({{out.outdir | quote}}, paste({{in.infile | fn | quote}}, "-", cname, {{in.infile | ext | quote}}, sep=""))
	write.table(m, ofile, sep = "\t", quote=F, row.names = rnames, col.names = cnames)
}