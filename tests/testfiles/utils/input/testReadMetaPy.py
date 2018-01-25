#!/usr/bin/env python

from bioprocs.utils import read

exec(read.meta.py)

rm = readMeta(('a', 'A'), ('b', 'B'))
assert isinstance(rm, readMeta)
assert rm.keys() == ['a', 'b']

rm.add(('c', 'C'))
rm.add(c = 'C')
assert rm.keys() == ['a', 'b', 'c']

rm.remove('b')
assert rm.keys() == ['a', 'c']

rm3 = readMeta()
rm3.borrow (rm)
assert rm3.keys() == ['a', 'c']

