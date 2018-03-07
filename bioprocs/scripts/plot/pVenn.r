mat         = read.table ({{in.infile | quote}}, sep="\t", header = TRUE, row.names = NULL, check.names = F)

{% if args.rnames %}
rns = make.unique(as.character(as.vector(mat[,1])))
mat[,1] = NULL
rownames(mat) = rns
{% endif %}

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