{{rimport}}('__init__.r')

infiles  = {{i.infiles | R}}
outfile  = {{o.outfile | R}}
params   = {{args.params | R}}
na       = {{args.na | R}}
fn2cname = {{args.fn2cname}}
fill     = as.logical({{args.fill | R}})
{% assign inopts = {} %}
{% for key, val in args.inopts.items() %}
	{% if isinstance(val, list) %}
		{% if len(val) == 1 %}
			{% python inopts[key] = val * len(i.infiles) %}
		{% else %}
			{% python inopts[key] = val %}
		{% endif %}
	{% else %}
		{% python inopts[key] = [val] * len(i.infiles) %}
	{% endif %}
{% endfor %}
inopts = {{inopts | R}}

mats = list()
for (i in 1:length(infiles)) {
	inparams = c(list(
		file      = infiles[i],
		header    = ifelse(is.null(inopts$cnames[i]), T, as.logical(inopts$cnames[i])),
		row.names = ifelse(is.null(inopts$rnames[i]) || inopts$rnames[i], 1, NULL),
		sep       = ifelse(is.null(inopts$delimit[i]), "\t", inopts$delimit[i]),
		skip      = ifelse(is.null(inopts$skip[i]), 0, inopts$skip[i])
	), params)
	mats[[i]]  = do.call(read.table, inparams)

	if (!inparams$header && !is.null(fn2cname)) {
		cname = fn2cname(tools::file_path_sans_ext(basename(infiles[i])))
		#if (ncol(mats[[i]]) == 1)
		#	colnames(mats[[i]]) = cname
		#else
		#	colnames(mats[[i]]) = paste(cname, 1:ncol(mats[[i]]), sep='_c')
		colnames(mats[[i]]) = cname
	}
}
mat = ifelse(fill, do.call(cbind.fill, mats), do.call(cbind, mats))
mat[is.na(mat)] = na

outparams = list(
	x         = mat,
	file      = outfile,
	sep       = '\t',
	quote     = F,
	col.names = T,
	row.names = ifelse(is.null(inopts$rnames[1]), T, as.logical(inopts$rnames[1]))
)
do.call(write.table, outparams)
