from pyppl import Box
from bioprocs.utils import regionOverlap
from bioprocs.utils.tsvio import TsvReader, TsvWriter

# read input file
regfile = {{i.regfile | quote}}
outfile = {{o.outfile | quote}}
ext     = {{i.regfile | ext | [1:] | quote}}
regopts = {{args.regopts}}
regions = None

if regopts.ftype == 'auto':
	regopts.ftype = ext
reader  = TsvReader(regfile, **regopts)
regions = reader.dump()

intfile = {{i.intfile | quote}}
ext     = {{i.intfile | ext | [1:] | quote}}
intopts = {{args.intopts}}

writer  = TsvWriter(outfile, ftype = 'bed')
writer.meta.add('CHR2', 'START2', 'END2', 'INTNAME')
writer.writeHead()

if intopts.ftype == 'auto':
	intopts.ftype = ext

import re
reader = TsvReader(intfile, **intopts)
if intopts.ftype == 'bedx':
	reader.meta.add('CHR2', 'START2', 'END2')

for line in reader:
	if intopts.ftype == 'bedpe':
		reg1 = [line.CHR1, line.START1, line.END1]
	else:
		reg1 = [line.CHR, line.START, line.END]
	if intopts.ftype == 'bed12':
		reg2 = re.split(r'[^\w]+', line.NAME)[:3]
	else: 
		reg2 = [line.CHR2, line.START2, line.END2]

	for region in regions:
		if regionOverlap(*(reg1 + [region.CHR, region.START, region.END])):
			region.CHR2   = reg2[0]
			region.START2 = reg2[1]
			region.END2   = reg2[2]
			region.INTNAME= line.NAME
			writer.write(region)
		if regionOverlap(*(reg2 + [region.CHR, region.START, region.END])):
			region.CHR2   = reg1[0]
			region.START2 = reg1[1]
			region.END2   = reg1[2]
			region.INTNAME= line.NAME
			writer.write(region)
