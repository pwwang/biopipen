from pyppl import Box
from os import path, symlink
from bioprocs.utils.gene import genenorm

genenorm(
	infile   = {{i.infile | quote}},
	outfile  = {{o.outfile | quote}},
	notfound = {{args.notfound | quote}},
	frm      = {{args.frm | quote}},
	to       = {{args.to  | quote}},
	genome   = {{args.genome | quote}},
	inopts   = {{args.inopts}},
	outopts  = {{args.outopts}},
	genecol  = {{args.genecol | quote}},
	cachedir = {{args.cachedir | quote}}
)
