from pyppl import Box
from bioprocs.utils.snp import snpinfo

infile   = {{in.snpfile | quote}}
outfile  = {{out.outfile | quote}}
genome   = {{args.genome | quote}}
dbsnpver = {{args.dbsnpver | quote}}
notfound = {{args.notfound | quote}}
inopts   = {{args.inopts}}
outopts  = {{args.outopts}}
snpcol   = {{args.snpcol | quote}}
cachedir = {{args.cachedir | quote}}

snpinfo (
	infile    = infile,
	outfile   = outfile,
	notfound  = notfound,
	genome    = genome,
	dbsnpver  = dbsnpver,
	inopts    = inopts,
	outopts   = outopts,
	snpcol    = snpcol,
	cachedir  = cachedir
)
