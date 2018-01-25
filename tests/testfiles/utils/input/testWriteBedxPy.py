#!/usr/bin/env python

from bioprocs.utils import write, read
from sys import argv
exec(read.bedx.py)
exec(write.bedx.py)

infile = argv[1]
outfile = argv[2]

reader = readBedx(infile)
writer = writeBedx(outfile)
writer.meta.borrow(reader.meta)

writer.writeMeta()
writer.writeHead()
for r in reader:
	writer.write(r)