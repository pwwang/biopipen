from pyppl import Box
from bioprocs.utils import shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile   = {{ i.infile | quote}}
outfile  = {{ o.outfile | quote}}
extend   = {{ args.extend | bool}}
gsize    = {{ args.gsize | quote}}
params   = {{ args.params | repr}}
bedtools = {{ args.bedtools | quote}}

shell.TOOLS.bedtools = bedtools
bedtools = shell.Shell(subcmd = True, dash = '-', equal = ' ').bedtools

params['g']   = gsize
params['i']   = infile

if not 'l' and not 'r' and not 'b' in params:
	raise ValueError('You have to define a length to flank (args.params.l, args.params.r or params.b')

if args.extend:
	left   = params.get('l', params.get('b', 0))
	right  = params.get('r', params.get('b', 0))
	stdns  = params.get('s', False)
	reader = TsvReader(infile, cnames = False)
	writer = TsvWriter(outfile)
	for r in reader:
		if not stdns or r[5] == '+':
			left2, right2 = left, right
		else:
			left2, right2 = right, left
		if params.pct:
			length   = r[2] - r[1]
			r[1] -= round(length * left2)
			r[2] += round(length * right2)
		else:
			r[1] -= left2
			r[2] += right2
		writer.write(r)
else:
	params._stdout = outfile
	bedtools.flank(**params).run()

