#!/usr/bin/env python

from bioprocs.utils import read
from sys import argv
exec(read.bed12.py)

infile = argv[1]

reader = readBed12(infile)
assert isinstance(reader, readBed12)

rs = reader.dump()
assert rs[0].CHR == 'chr22'
assert rs[0].START == '1000'
assert rs[0].END == '5000'
assert rs[0].NAME == 'cloneA'
assert rs[0].SCORE == '960'
assert rs[0].STRAND == '+'
assert rs[0].THICKSTART == '1000'
assert rs[0].THICKEND == '5000'
assert rs[0].ITEMRGB == '0'
assert rs[0].BLOCKCOUNT == '2'
assert rs[0].BLOCKSIZES == '567,488,'
assert rs[0].BLOCKSTARTS == '0,3512'
assert rs[1].CHR == 'chr22'
assert rs[1].START == '2000'
assert rs[1].END == '6000'
assert rs[1].NAME == 'cloneB'
assert rs[1].SCORE == '900'
assert rs[1].STRAND == '-'
assert rs[1].THICKSTART == '2000'
assert rs[1].THICKEND == '6000'
assert rs[1].ITEMRGB == '0'
assert rs[1].BLOCKCOUNT == '2'
assert rs[1].BLOCKSIZES == '433,399,'
assert rs[1].BLOCKSTARTS == '0,3601'


