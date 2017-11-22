cnames = as.logical({{args.cnames | R}})
rnames = as.logical({{args.rnames | R}})

mat    = read.table ({{in.infile | quote}}, sep="\t", header = cnames, row.names = if (rnames) 1 else NULL, check.names = F)
rns    = if (rnames) rownames(mat) else paste("ROW", 1:nrow(mat), sep="")
for (i in 1:nrow(mat)) {
	rname = rns[i]
	m     = mat[rname, , drop=F]
	ofile = file.path({{out.outdir | quote}}, paste({{in.infile | fn | quote}}, "-", rname, {{in.infile | ext | quote}}, sep=""))
	write.table(m, ofile, sep = "\t", quote=F, row.names = rnames, col.names = cnames)
}