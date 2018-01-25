#!/usr/bin/env python

from bioprocs.utils import read

exec(read.record.py)

rr = readRecord(a='A', b='B')

assert isinstance(rr, readRecord)
assert sorted(rr.keys()) == ['a', 'b']

rr.add(c = 'C')
assert sorted(rr.keys()) == ['a', 'b', 'c']

