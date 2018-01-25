#!/usr/bin/env python

from bioprocs.utils import read
from sys import argv
exec(read.base.py)

infile = argv[1]

reader = readBase(infile)
reader.meta = readMeta("A", "B", "C", "D", "E")
assert isinstance(reader, readBase)

for r in reader:
	assert r.A == 'a1'
	assert r.B == 'b1'
	assert r.C == 'c1'
	assert r.D == 'd1'
	assert r.E == 'e1'
	break

reader.rewind()
rs = reader.dump()
assert rs[0].A == 'a1'
assert rs[0].B == 'b1'
assert rs[0].C == 'c1'
assert rs[0].D == 'd1'
assert rs[0].E == 'e1'
assert rs[1].A == 'a2'
assert rs[1].B == 'b2'
assert rs[1].C == 'c2'
assert rs[1].D == 'd2'
assert rs[1].E == 'e2'


