cnames = as.logical({{args.cnames | R}})
rnames = as.logical({{args.rnames | R}})

mat    = read.table ({{in.infile | quote}}, sep="\t", header = cnames, row.names = if (rnames) 1 else NULL, check.names = F)
if (cnames) {
	heads = colnames(mat)
} else {
	heads = paste("COL", 1:ncol(mat), sep="")
	colnames(mat) = heads
}
n = {{args.n}}
chunkfunc = function(d) split(d, ceiling(seq_along(d)/n))
chunks = chunkfunc(1:ncol(mat))
for (chunk in chunks) {
	cname = heads[chunk]
	m     = mat[, cname, drop=F]
	ofile = file.path({{out.outdir | quote}}, paste({{in.infile | fn | quote}}, "-", gsub('[^a-zA-Z0-9]', '', cname[1]), {{in.infile | ext | quote}}, sep=""))
	write.table(m, ofile, sep = "\t", quote=F, row.names = rnames, col.names = cnames)
}