{{plotPie}}

ggs = {{args.ggs | Rlist}}
mat = read.table({{in.infile | quote}}, sep = "\t", header = T, row.names = NULL)
{% if args.rnames %}
rns = make.unique(as.character(as.vector(mat[,1])))
mat[,1] = NULL
rownames(mat) = rns
{% endif %}
plotPie(mat, filename = {{out.outfile | quote}}, ggs = ggs, devpars = {{args.devpars | Rlist}})