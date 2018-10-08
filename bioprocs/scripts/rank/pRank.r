{{rimport}}('__init__.r')

infile   = {{i.infile    | quote}}
outfile  = {{o.outfile  | quote}}
byrow    = {{args.byrow   | R}}
inopts   = {{args.inopts  | R}}
method   = {{args.tie     | quote}}
na       = {{args.na      | quote}}
reverse  = {{args.reverse | R}}

inparams = list(
	file        = infile,
	header      = ifelse(is.null(inopts$cnames) || inopts$cnames, T, F),
	row.names   = ifelse(is.null(inopts$rnames) || inopts$rnames, 1, NULL),
	sep         = ifelse(is.null(inopts$delimit), "\t", inopts$delimit),
	check.names = F
)

data = do.call(read.table.nodup, inparams)
if (reverse) data = -data
if (byrow)   data = t(data)

if (na == 'remove') {
	na = NA
} else if (na == 'first') {
	na = F
} else if (na == 'last') {
	na = T
} # else keep

ranks = apply(data, 2, function(row) rank(row, na.last = na, ties.method = method))
if (byrow) ranks = t(ranks)

outparams = list(
	x         = ranks,
	file      = outfile,
	quote     = F,
	sep       = inparams$sep,
	row.names = ifelse(byrow, is.null(inopts$cnames) || inopts$cnames, is.null(inopts$rnames) || inopts$rnames),
	col.names = ifelse(byrow, is.null(inopts$rnames) || inopts$rnames, is.null(inopts$cnames) || inopts$cnames)
)
do.call(write.table, outparams)