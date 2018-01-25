#!/usr/bin/env python

from bioprocs.utils import read
from sys import argv
exec(read.bedpe.py)

infile = argv[1]

reader = readBedpe(infile)
assert isinstance(reader, readBedpe)

rs = reader.dump()
assert rs[0].CHR1 == 'chr1'
assert rs[0].START1 == '100'
assert rs[0].END1 == '200'
assert rs[0].CHR2 == 'chr5'
assert rs[0].START2 == '5000'
assert rs[0].END2 == '5100'
assert rs[0].NAME == 'bedpe_example1'
assert rs[0].SCORE == '30'
assert rs[0].STRAND1 == '+'
assert rs[0].STRAND2 == '-'
assert rs[1].CHR1 == 'chr9'
assert rs[1].START1 == '1000'
assert rs[1].END1 == '5000'
assert rs[1].CHR2 == 'chr9'
assert rs[1].START2 == '3000'
assert rs[1].END2 == '3800'
assert rs[1].NAME == 'bedpe_example2'
assert rs[1].SCORE == '100'
assert rs[1].STRAND1 == '+'
assert rs[1].STRAND2 == '-'