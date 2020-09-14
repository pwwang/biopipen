"""Script for bedtools.pBedFlank"""
# pylint: disable=invalid-name,undefined-variable,unused-import
# pylint: disable=unsupported-assignment-operation,not-a-mapping
from diot import Diot
from bioprocs.utils import shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
extend = {{args.extend | bool}}
gsize = {{args.gsize | quote}}
params = {{args.params | repr}}
bedtools = {{args.bedtools | quote}}

shell.load_config(bedtools=bedtools)

params['g'] = gsize
params['i'] = infile

if not 'l' and not 'r' and not 'b' in params:
    raise ValueError('You have to define a length to flank '
                     '(args.params.l, args.params.r or params.b')

if extend:
    left = params.get('l', params.get('b', 0))
    right = params.get('r', params.get('b', 0))
    stdns = params.get('s', False)
    reader = TsvReader(infile, cnames=False)
    writer = TsvWriter(outfile)
    for r in reader:
        r[1] = int(r[1])
        r[2] = int(r[2])
        if not stdns or r[5] == '+':
            left2, right2 = left, right
        else:
            left2, right2 = right, left
        if params.get('pct'):
            length = r[2] - r[1]
            r[1] -= round(length * left2)
            r[2] += round(length * right2)
        else:
            r[1] -= left2
            r[2] += right2
        r[1] = max(0, r[1])
        writer.write(r)
else:
    # params._out = outfile
    # params._debug = True
    shell.bedtools.flank(**params).r > outfile
