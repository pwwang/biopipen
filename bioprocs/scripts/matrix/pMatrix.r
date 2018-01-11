cnames = as.logical({{args.cnames | R}})
rnames = as.logical({{args.rnames | R}})

infile           = {{in.infile | quote}}
params           = {{args.params | Rlist}}
params$file      = infile
params$header    = cnames
params$row.names = if (rnames) 1 else NULL
mat              = do.call(read.table, params)
if (rnames) {
	rns = make.unique(as.vector(mat[,1]))
	mat[,1] = NULL
	rownames(mat) = rns
}
{% if args.code | lambda x: isinstance(x, list) %}
{% for c in args.code %}
{{c}}
{% endfor %}
{% else %}
{{args.code}}
{% endif %}
write.table (mat, {{out.outfile | quote}}, sep=params$sep, quote=F, col.names = cnames, row.names = rnames)
