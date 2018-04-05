from bioprocs.utils import regionOverlap
from bioprocs.utils.tsvio import TsvReader, TsvWriter

# read input file
regfile = {{in.regfile | quote}}
outfile = {{out.outfile | quote}}
ext     = {{in.regfile | ext | [1:] | quote}}
regtype = {{args.regtype | quote}}
regions = None

if (regtype == 'auto' and ext == 'bed') or regtype == 'bed':
	reader  = TsvReader(regfile, ftype = 'bed')
	regions = reader.dump()
elif (regtype == 'auto' and ext == 'bedx') or regtype == 'bedx':
	reader  = TsvReader(regfile, ftype = 'bedx')
	regions = reader.dump()

intfile = {{in.intfile | quote}}
ext     = {{in.intfile | ext | [1:] | quote}}
inttype = {{args.inttype | quote}}

writer  = TsvWriter(outfile, ftype = 'bed')
writer.meta.add('CHR2', 'START2', 'END2', 'INTNAME')
writer.writeHead()
if (inttype == 'auto' and ext == 'bed12') or inttype == 'bed12':
	import re
	reader = TsvReader(intfile, ftype = 'bed12')
	for line in reader:
		reg1 = [line.CHR, line.START, line.END]
		reg2 = re.split(r'[^\w]+', line.NAME)[:3]
		for region in regions:
			if regionOverlap(*(reg1 + [region.CHR, region.START, region.END])):
				region.CHR2   = reg2[0]
				region.START2 = reg2[1]
				region.END2   = reg2[2]
				region.INTNAME= line.NAME
				writer.write(region)
			if regionOverlap(*(reg2 + [region.CHR, region.START, region.END])):
				region.CHR2   = line.CHR
				region.START2 = line.START
				region.END2   = line.END
				region.INTNAME= line.NAME
				writer.write(region)
elif (inttype == 'auto' and ext == 'bedx') or inttype == 'bedx':
	reader = TsvReader(intfile, ftype = 'bedx')
	reader.meta.add('CHR2', 'START2', 'END2')
	for line in reader:
		reg1 = [line.CHR, line.START, line.END]
		reg2 = [line.CHR2, line.START2, line.END2]
		for region in regions:
			if regionOverlap(*(reg1 + [region.CHR, region.START, region.END])):
				region.CHR2   = line.CHR2
				region.START2 = line.START2
				region.END2   = line.END2
				region.INTNAME= line.NAME
				writer.write(region)
			if regionOverlap(*(reg2 + [region.CHR, region.START, region.END])):
				region.CHR2   = line.CHR
				region.START2 = line.START
				region.END2   = line.END
				region.INTNAME= line.NAME
				writer.write(region)

elif (inttype == 'auto' and ext == 'bedpe') or inttype == 'bedpe':
	reader = TsvReader(intfile, ftype = 'bedpe')
	for line in reader:
		reg1 = [line.CHR1, line.START1, line.END1]
		reg2 = [line.CHR2, line.START2, line.END2]
		for region in regions:
			if regionOverlap(*(reg1 + [region.CHR, region.START, region.END])):
				region.CHR2   = line.CHR2
				region.START2 = line.START2
				region.END2   = line.END2
				region.INTNAME= line.NAME
				writer.write(region)
			if regionOverlap(*(reg2 + [region.CHR, region.START, region.END])):
				region.CHR2   = line.CHR1
				region.START2 = line.START1
				region.END2   = line.END1
				region.INTNAME= line.NAME
				writer.write(region)
			
			
			





