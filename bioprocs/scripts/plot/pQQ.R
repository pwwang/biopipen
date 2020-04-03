{{"__init__.R", "plot.R" | rimport}}
options(stringsAsFactors = FALSE)

infile  = {{ i.infile | quote}}
outfile = {{ o.outfile | quote}}
inopts  = {{ args.inopts | R}}
devpars = {{ args.devpars | R}}
ggs     = {{ args.ggs | R}}
tsform  = {{ args.tsform | ?.startswith: "function" | !R }}
params  = {{ args.params | R}}

indata = read.table.inopts(infile, inopts)

plot.qq(
	indata,
	plotfile = outfile,
	x        = ifelse(ncol(indata) == 1, NULL, 1),
	y        = ifelse(ncol(indata) == 1, 1, 2),
	params   = c(list(tsform = tsform), params),
	ggs      = ggs,
	devpars  = devpars
)
