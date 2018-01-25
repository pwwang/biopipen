#!/usr/bin/env python

from bioprocs.utils import write, read
from sys import argv
exec(read.base.py)
exec(write.base.py)

infile = argv[1]
outfile = argv[2]

reader = readBase(infile)
reader.meta = readMeta("A", "B", "C", "D", "E")
writer = writeBase(outfile)
writer.meta.borrow(reader.meta)
writer.meta.remove('C')

writer.writeMeta()
writer.writeHead()
for r in reader:
	writer.write(r)