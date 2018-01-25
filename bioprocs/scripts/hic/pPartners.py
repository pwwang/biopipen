def overlap(CHR1, START1, END1, CHR2, START2, END2):
	if CHR1 != CHR2: return False
	if int(END1) < int(START2): return False
	if int(END2) < int(START1): return False
	return True

{{readBedx}}
{{readBed}}
{{readBed12}}
{{readBedpe}}
{{writeBedx}}
# read input file
ext     = {{in.regfile | ext | [1:] | quote}}
regtype = {{args.regtype | quote}}
regions = None

if (regtype == 'auto' and ext == 'bed') or regtype == 'bed':
	reader = readBed({{in.regfile | quote}})
	regions = reader.dump()
elif (regtype == 'auto' and ext == 'bedx') or regtype == 'bedx':
	reader = readBedx({{in.regfile | quote}})
	regions = reader.dump()

ext     = {{in.intfile | ext | [1:] | quote}}
inttype = {{args.inttype | quote}}

writer  = writeBedx({{out.outfile | quote}})
writer.meta.add(*readBedx.META)
writer.meta.add(
	('CHR2', 'Interaction partner chromosome'),
	('START2', 'Interaction partner start'),
	('END2', 'Interaction partner end'),
	('INTNAME', 'Interaction name')
)
writer.writeMeta()
writer.writeHead()
if (inttype == 'auto' and ext == 'bed12') or inttype == 'bed12':
	import re
	reader = readBed({{in.intfile | quote}})
	for line in reader:
		reg1 = [line.CHR, line.START, line.END]
		reg2 = re.split(r'[^\w]+', line.NAME)[:3]
		for region in regions:
			if overlap(*(reg1 + [region.CHR, region.START, region.END])):
				region.CHR2   = reg2[0]
				region.START2 = reg2[1]
				region.END2   = reg2[2]
				region.INTNAME= line.NAME
				writer.write(region)
			if overlap(*(reg2 + [region.CHR, region.START, region.END])):
				region.CHR2   = line.CHR
				region.START2 = line.START
				region.END2   = line.END
				region.INTNAME= line.NAME
				writer.write(region)
elif (inttype == 'auto' and ext == 'bedx') or inttype == 'bedx':
	reader = readBedx({{in.intfile | quote}})
	for line in reader:
		reg1 = [line.CHR, line.START, line.END]
		reg2 = [line.CHR2, line.START2, line.END2]
		for region in regions:
			if overlap(*(reg1 + [region.CHR, region.START, region.END])):
				region.CHR2   = line.CHR2
				region.START2 = line.START2
				region.END2   = line.END2
				region.INTNAME= line.NAME
				writer.write(region)
			if overlap(*(reg2 + [region.CHR, region.START, region.END])):
				region.CHR2   = line.CHR
				region.START2 = line.START
				region.END2   = line.END
				region.INTNAME= line.NAME
				writer.write(region)

elif (inttype == 'auto' and ext == 'bedpe') or inttype == 'bedpe':
	reader = readBedpe({{in.intfile | quote}})
	for line in reader:
		reg1 = [line.CHR1, line.START1, line.END1]
		reg2 = [line.CHR2, line.START2, line.END2]
		for region in regions:
			if overlap(*(reg1 + [region.CHR, region.START, region.END])):
				region.CHR2   = line.CHR2
				region.START2 = line.START2
				region.END2   = line.END2
				region.INTNAME= line.NAME
				writer.write(region)
			if overlap(*(reg2 + [region.CHR, region.START, region.END])):
				region.CHR2   = line.CHR1
				region.START2 = line.START1
				region.END2   = line.END1
				region.INTNAME= line.NAME
				writer.write(region)
			
			
			





