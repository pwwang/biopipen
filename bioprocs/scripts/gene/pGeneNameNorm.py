from os import path, symlink
{{ genenorm | norepeats }}

_, cachefile = genenorm(
	{{in.infile | quote}}, 
	outfile  = {{out.outfile | quote}}, 
	notfound = {{args.notfound | quote}},
	frm      = {{args.frm | quote}},
	to       = {{args.to  | quote}},
	genome   = {{args.genome | quote}},
	tmpdir   = {{args.tmpdir | quote}},
	inopts   = {{args.inopts}},
	outopts  = {{args.outopts}},
	inmeta   = {{args.inmeta | lambda x: x if isinstance(x, list) or isinstance(x, dict) else '"' + x + '"'}},
	genecol  = {{args.genecol | quote}}
)

symlink(cachefile, path.join({{job.outdir | quote}}, path.basename(cachefile)))