mat         = read.table ({{in.infile | quote}}, sep="\t", header = TRUE, row.names = {{args.rnames | lambda x: '1' if x else 'NULL'}}, check.names = F)
nc          = ncol(mat)
cnames      = colnames(mat)
vennParams  = {{args.vennParams  | Rlist}}
upsetParams = {{args.upsetParams | Rlist}}
devpars     = {{args.devpars     | Rlist}}
# use real venn plot 
if (nc<=3 && ({{args.tool | quote}} == 'auto' || {{args.tool | quote}} == 'venn')) {
	{{plotVenn}}
	plotVenn (mat, {{out.outfile | quote}}, vennParams, devpars)

# use upset plot #
} else {
	{{plotUpset}}
	plotUpset(mat, {{out.outfile | quote}}, upsetParams, devpars)
}