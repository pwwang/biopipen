from pyppl import Box
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile1  = {{i.infile1 | quote}}
infile2  = {{i.infile2 | quote}}
outfile1 = {{o.outfile1 | quote}}
outfile2 = {{o.outfile2 | quote}}
inopts1  = {{args.inopts1 | repr}}
inopts2  = {{args.inopts2 | repr}}

if inopts1.headCallback:
	inopts1.headCallback = eval(inopts1.headCallback)
if inopts2.headCallback:
	inopts2.headCallback = eval(inopts2.headCallback)
inopts1.attach = False
inopts2.attach = False

rnames1 = inopts1.get('rnames', True)
rnames2 = inopts2.get('rnames', True)
if 'rnames' in inopts1:
	del inopts1['rnames']
if 'rnames' in inopts2:
	del inopts2['rnames']

indata1 = TsvReader(infile1, **inopts1)
indata2 = TsvReader(infile2, **inopts2)

cnames1 = indata1.meta if not rnames1 else indata1.meta[1:]
cnames2 = indata2.meta if not rnames2 else indata2.meta[1:]
paired  = list(set(cnames1) & set(cnames2))
cnames1 = cnames2 = paired

if rnames1:
	cnames1 = [indata1.meta[0]] + cnames1
if rnames2:
	cnames2 = [indata2.meta[0]] + cnames2

cindex1 = [indata1.meta.index(c) for c in cnames1]
cindex2 = [indata2.meta.index(c) for c in cnames2]

outdata1 = TsvWriter(outfile1)
outdata2 = TsvWriter(outfile2)
outdata1.meta = cnames1
outdata2.meta = cnames2
outdata1.writeHead()
outdata2.writeHead()

for r1 in indata1:
	outdata1.write(r1[i] for i in cindex1)
outdata1.close()

for r2 in indata2:
	outdata2.write(r2[i] for i in cindex2)
outdata2.close()
