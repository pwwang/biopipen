{{ '__init__.R' | rimport }}

library(tidyr)
library(dplyr)

infile  = {{in.infile | r}}
outfile = {{out.outfile | r}}
inopts  = {{args.inopts | r}}

df = read.table.inopts(infile, inopts)

{{args.code | isinstance: list ? join: '\n' ! | render }}

outparams = list(
	x         = df,
	file      = outfile,
	sep       = list.get(inopts, "delimit", "\t"),
	quote     = F,
	col.names = is.true(list.get(inopts, 'cnames', TRUE)),
	row.names = is.true(list.get(inopts, 'rnames', TRUE))
)

do.call(write.table, outparams)
