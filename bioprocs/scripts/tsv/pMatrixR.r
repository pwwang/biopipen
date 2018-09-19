{{rimport}}('__init__.r')

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
params  = {{args.params | R}}
inopts  = {{args.inopts | R}}

inparams = list(
	file      = infile,
	header    = as.logical(inopts$cnames),
	row.names = if (as.logical(inopts$rnames)) 1 else NULL,
	skip      = if (is.null(inopts$skip)) 0 else as.numeric(inopts$skip),
	sep       = inopts$delimit
)

mat = do.call(read.table.nodup, c(inparams, params))

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
