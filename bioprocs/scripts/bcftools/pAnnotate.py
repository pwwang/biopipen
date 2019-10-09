
from os import path
from pyppl import Box
from bioprocs.utils import shell2 as shell
from bioprocs.utils.reference import tabixIndex
from tempfile import gettempdir

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params | repr}}
nthread  = {{args.nthread | repr}}
bcftools = {{args.bcftools | quote}}
annfile  = {{args.annfile | quote}}
tabix    = {{args.tabix | quote}}
header   = {{args.header | repr}}
cols     = {{args.cols | repr}}

shell.load_config(bcftools = bcftools)

if annfile and path.isfile(annfile):
	abname   = path.basename(annfile)
	ext      = path.splitext(abname[:-3] if abname.endswith('.gz') else abname)[-1][1:]
	infile   = tabixIndex(infile, 'vcf', tabix)
	params.a = tabixIndex(annfile, ext, tabix)
if cols:
	params.c = ','.join(cols if isinstance(cols, list) else [c for c in cols.split(',')])
if header:
	if not isinstance(header, list):
		header = [header]
	from datetime import datetime
	now = str(datetime.now().timestamp())
	headerfile = path.join(gettempdir(), 'pBcftools.annotate.header.%s' % now)
	with open(headerfile, 'w') as fhead:
		for head in header:
			fhead.write(head + '\n')
	params.h = headerfile

params._ = infile 
params.o = outfile
params.threads = nthread

shell.fg.bcftools.annotate(**params)
