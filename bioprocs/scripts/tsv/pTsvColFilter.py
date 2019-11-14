from os import path
from collections import OrderedDict
from pyppl import Box
from pyppl.utils import alwaysList
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile  = {{i.infile | quote}}
colfile = {{i.colfile | quote}}
outfile = {{o.outfile | quote}}
inopts  = {{args.inopts | repr}}
cols    = {{args.cols | repr}}
keep    = {{args.keep | repr}}

if path.isfile(colfile):
	cols = TsvReader(colfile, cnames = False).dump(0)
elif colfile:
	cols = alwaysList(colfile)
elif path.isfile(str(cols)):
	cols = TsvReader(cols, cnames = False).dump(0)
elif cols:
	cols = alwaysList(cols)
else:
	raise ValueError('Columns not provided.')
if not isinstance(cols[0], int) and cols[0].isdigit():
	cols = [int(c) for c in cols]

reader = TsvReader(infile, **inopts)
writer = TsvWriter(outfile, delimit = inopts.get('delimit', "\t"))
if reader.cnames and not isinstance(cols[0], int):
	cols = [reader.cnames.index(c) for c in cols if c in reader.cnames]
elif not reader.cnames and not isinstance(cols[0], int):
	raise ValueError("Input file doesn't have column names")
elif min(cols) < -len(reader.cnames) or (reader.cnames and max(cols) >= len(reader.cnames)):
	raise IndexError("Provided columns beyond input file range.")

if reader.cnames:
	ncol = len(reader.cnames)
else:
	ncol = len(reader.next())
	reader.rewind()

cols = [ncol + c if c < 0 else c for c in cols]
if not keep:
	cols = [c for c in range(ncol) if c not in cols]

if reader.cnames:
	writer.cnames = [reader.cnames[c] for c in cols]
	writer.writeHead()

for r in reader:
	rec = [r[c] for c in cols]
	writer.write(rec)
writer.close()
