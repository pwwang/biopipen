from pyppl import Box
from bioprocs.utils.helpers import runcmd, cmdargs

params = Box()
params['g'] = {{args.gsize | quote}}
params['i'] = {{in.infile | quote}}
params['pct'] = False
params.update({{args.params}})

if not 'l' and not 'r' and not 'b' in params:
	raise ValueError('You have to define a length to flank (args.params.l, args.params.r or params.b')

{% if args.extend %}
left  = params['l'] if 'l' in params else params['b'] if 'b' in params else 0
right = params['r'] if 'r' in params else params['b'] if 'b' in params else 0
stdns = params['s'] if 's' in params else False # strandness
from bioprocs.utils.tsvio import TsvReader, TsvWriter
reader = TsvReader({{in.infile | quote}}, ftype = 'bed')
writer = TsvWriter({{out.outfile | quote}})
writer.meta.update(reader.meta)
for r in reader:
	if not stdns or r.STRAND == '+':
		left2, right2 = left, right
	else:
		left2, right2 = right, left
	if params['pct']:
		length   = r.END - r.START
		r.START -= round(length * left2)
		r.END   += round(length * right2)
	else:
		r.START -= left2
		r.END   += right2
	writer.write(r)
{% else %}
cmd = '{{args.bedtools}} flank %s > {{out.outfile | quote}}' % cmdargs(params, dash='-', equal=' ')
runcmd(cmd)
{% endif %}