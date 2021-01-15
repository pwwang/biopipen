{{'__init__.R' | rimport}}

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
params  = {{args.params | R}}
inopts  = {{args.inopts | R}}

mat = read.table.inopts(infile, inopts)

{{args.code | ?isinstance: list | = '\n'.join | $ render }}

outparams = list(
	x         = mat,
	file      = outfile,
	sep       = inopts$delimit,
	quote     = F,
	col.names = as.logical(inopts$cnames),
	row.names = as.logical(inopts$rnames)
)

do.call(write.table, outparams)
