#!/usr/bin/env python

from bioprocs.utils import read
from sys import argv
exec(read.bedx.py)

infile = argv[1]

reader = readBedx(infile)
assert isinstance(reader, readBedx)

rs = reader.dump()
assert rs[0].CHR == 'chr1'
assert rs[0].START == '100'
assert rs[0].END == '200'
assert rs[0].CHR2 == 'chr5'
assert rs[0].START2 == '5000'
assert rs[0].END2 == '5100'
assert rs[0].NAME == 'bedx1'
assert rs[0].SCORE == '30'
assert rs[0].STRAND == '+'
assert rs[0].STRAND2 == '-'
assert rs[1].CHR == 'chr9'
assert rs[1].START == '1000'
assert rs[1].END == '5000'
assert rs[1].CHR2 == 'chr9'
assert rs[1].START2 == '3000'
assert rs[1].END2 == '3800'
assert rs[1].NAME == 'bedx2'
assert rs[1].SCORE == '100'
assert rs[1].STRAND == '+'
assert rs[1].STRAND2 == '-'