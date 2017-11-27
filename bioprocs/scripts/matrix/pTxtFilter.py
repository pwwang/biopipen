{{txtFilter}}

txtFilter(
	infile  = {{in.infile | quote}},
	outfile = {{out.outfile | quote}},
	cols    = {{args.cols}},
	rfilter = {{args.rfilter}},
	header  = {{ args.header }},
	skip    = {{args.skip}},
	delimit = {{args.delimit | quote}},
	outdelimit = {{args.outdelimit | quote}}
)