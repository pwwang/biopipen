{{plotHeatmap}}

m = read.table({{in.infile | quote}}, sep="\t", header=T, row.names=1, check.names=F)
{% if args | lambda x: not 'dendro' in x %}
dendro = T
rows   = NULL
cols   = NULL
{% else %}
dendro = {{args.dendro.dendro | R}}
rows   = {% if args.dendro | lambda x: 'rows' in x %}{{args.dendro.rows | Rvec}}{% else %}NULL{% endif %}
cols   = {% if args.dendro | lambda x: 'cols' in x %}{{args.dendro.cols | Rvec}}{% else %}NULL{% endif %}
{% endif %}
plotHeatmap(m, {{out.outfile | quote}}, ggs={{args | lambda x: {} if not 'ggs' in x else x['ggs'] | Rlist}}, dendro = dendro, rows = rows, cols = cols)
