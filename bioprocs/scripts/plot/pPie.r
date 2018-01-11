{{plotPie}}

ggs = {{args.ggs | Rlist}}
mat = read.table({{in.infile | quote}}, sep = "\t", header = T, row.names = {% if args.rnames %}1{% else %}NULL{% endif %})
plotPie(mat, filename = {{out.outfile | quote}}, ggs = ggs, devpars = {{args.devpars | Rlist}})