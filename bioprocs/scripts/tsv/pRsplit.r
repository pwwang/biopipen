cnames = as.logical({{args.cnames | R}})
rnames = as.logical({{args.rnames | R}})

mat    = read.table ({{in.infile | quote}}, sep="\t", header = cnames, row.names = if (rnames) 1 else NULL, check.names = F)
if (rnames) {
	rns = rownames(mat)
} else {
	rns = paste("ROW", 1:nrow(mat), sep="")
	rownames(mat) = rns
}
n = {{args.n}}
chunkfunc = function(d) split(d, ceiling(seq_along(d)/n))
chunks = chunkfunc(1:nrow(mat))
for (chunk in chunks) {
	rname = rns[chunk]
	m     = mat[rname, , drop=F]
	ofile = file.path({{out.outdir | quote}}, paste({{in.infile | fn | quote}}, "-", gsub('[^a-zA-Z0-9]', '', rname[1]), {{in.infile | ext | quote}}, sep=""))
	write.table(m, ofile, sep = "\t", quote=F, row.names = rnames, col.names = cnames)
}