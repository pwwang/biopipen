cnames = as.logical({{args.cnames | R}})
rnames = as.logical({{args.rnames | R}})

infile           = {{in.infile | quote}}
params           = {{args.params | Rlist}}
params$file      = infile
params$header    = cnames
mat              = do.call(read.table, c(list(row.names = NULL), params))
if (rnames) {
	rns = as.vector(mat[,1])
	if (length(rns) > 0) {
		rns = make.unique(rns)
		mat[,1] = NULL
		rownames(mat) = rns
	}
}
{% if args.code | lambda x: isinstance(x, list) %}
{% for c in args.code %}
{{c}}
{% endfor %}
{% else %}
{{args.code}}
{% endif %}
write.table (mat, {{out.outfile | quote}}, sep=params$sep, quote=F, col.names = cnames, row.names = rnames)
