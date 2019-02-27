{{rimport}}('__init__.r')

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
params  = {{args.params | R}}
inopts  = {{args.inopts | R}}

mat = read.table.inopts(infile, inopts)

{% if isinstance(args.code, list) %}
	{% for c in args.code %}
	{{c}}
	{% endfor %}
{% else %}
	{{args.code}}
{% endif %}

outparams = list(
	x         = mat,
	file      = outfile,
	sep       = inopts$delimit,
	quote     = F,
	col.names = as.logical(inopts$cnames),
	row.names = as.logical(inopts$rnames)
)

do.call(write.table, outparams)
