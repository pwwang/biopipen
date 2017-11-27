{{txtTransform}}

txtTransform(
	infile    = {{in.infile | quote}},
	outfile   = {{out.outfile | quote}},
	cols      = {{args.cols}},
	transform = {{args.transform}},
	header    = {{ args.header }},
	skip      = {{args.skip}},
	delimit   = {{args.delimit | quote}},
	outdelimit= {{args.outdelimit | quote}}
)