mat      = read.table ("{{in.infile}}", sep="\t", header = {{args.header | R}}, row.names = {{args.rownames | R }}, check.names = F)
{% if args.code | lambda x: isinstance(x, list) %}
{% for c in args.code %}
{{c}}
{% endfor %}
{% else %}
{{args.code}}
{% endif %}
write.table (mat, {{out.outfile | quote}}, sep="\t", quote=F, col.names = {{args.header | R}}, row.names = as.logical({{args.rownames | R}}))
