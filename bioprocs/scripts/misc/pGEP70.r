{{rimport}}('__init__.r')
library(methods)

infile  = {{in.infile | R}}
outfile = {{out.outfile | R}}
gep70   = {{args.gep70 | R}}
inopts  = {{args.inopts | R}}

inopts.default = list(
	cnames  = T,
	rnames  = T,
	delimit = "\t",
	skip    = 0
)
inopts = update.list(inopts.default, inopts)

# read all data
inparams = list(
	sep       = inopts$delimit,
	file      = infile,
	header    = inopts$cnames,
	row.names = ifelse(inopts$rnames, 1, NULL),
	skip      = inopts$skip
)
mat = do.call(read.table, inparams)

# read the 70 genes
gene70 = read.table(gep70, row.names = NULL, header = F, check.names = F, fill = T)

allgenes = rownames(mat)
gene51   = intersect(allgenes, gene70[,1])
gene19   = intersect(allgenes, gene70[,2])
exp51    = colMeans(mat[gene51, , drop = F])
exp19    = colMeans(mat[gene19, , drop = F])
ret      = exp51 - exp19

rownames(ret) = 'GEP70'
ret = t(ret)
write.table(ret, outfile, quote = F, row.names = inopts$cnames, col.names = T, sep = "\t")


