cnames = as.logical({{args.cnames | R}})
rnames = as.logical({{args.rnames | R}})

mat      = read.table ("{{in.infile}}", sep="\t", header = cnames, row.names = if (rnames) 1 else NULL, check.names = F)
{% if args.code | lambda x: isinstance(x, list) %}
{% for c in args.code %}
{{c}}
{% endfor %}
{% else %}
{{args.code}}
{% endif %}
write.table (mat, {{out.outfile | quote}}, sep="\t", quote=F, col.names = cnames, row.names = rnames)
