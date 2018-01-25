#!/usr/bin/env python

from bioprocs.utils import read
from sys import argv
exec(read.bed.py)

infile = argv[1]

reader = readBed(infile)
assert isinstance(reader, readBed)

rs = reader.dump()
assert rs[0].CHR == 'chr1'
assert rs[0].START == '100'
assert rs[0].END == '200'
assert rs[0].NAME == 'bedx1'
assert rs[0].SCORE == '30'
assert rs[0].STRAND == '+'
assert rs[1].CHR == 'chr9'
assert rs[1].START == '1000'
assert rs[1].END == '5000'
assert rs[1].NAME == 'bedx2'
assert rs[1].SCORE == '100'
assert rs[1].STRAND == '+'