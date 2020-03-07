"""Script for bed.pBedSwitchBase"""

from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
inbase = {{arg.inbase | repr}}
outbase = {{arg.outbase | repr}}

if inbase is None and outbase is None:
    raise ValueError("args.inbase or args.outbase has to be specified.")

if inbase is None:
    inbase = 1 - outbase

if outbase is None:
    outbase = 1 - inbase

if inbase not in (0, 1) or outbase not in (0, 1):
    raise ValueError("args.inbase and args.outbase have to be 0 or 1.")

if inbase == outbase:
    raise ValueError("Not transforming to the same base")

def row_transform(row):
    """Transform a BED row between bases"""
    start = int(row[1])
    if inbase == 0:
        start += 1
    else:
        start -= 1
    row[1] = start
    return row

reader = TsvReader(infile, cnames=False)
writer = TsvWriter(outfile)

for row in reader:
    writer.write(row_transform(row))

reader.close()
writer.close()
