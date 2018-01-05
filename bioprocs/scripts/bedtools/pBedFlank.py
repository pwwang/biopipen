{{runcmd}}
{{params2CmdArgs}}

params = {{args.params}}
params['g'] = {{args.gsize | quote}}
params['i'] = {{in.infile | quote}}

if not 'l' and not 'r' and not 'b' in params:
	raise ValueError('You have to define a length to flank (args.params.l, args.params.r or params.b')

{% if args.extend %}
left  = params['l'] if 'l' in params else params['b'] if 'b' in params else 0
right = params['r'] if 'r' in params else params['b'] if 'b' in params else 0
stdns = params['s'] if 's' in params else False # strandness
with open({{in.infile | quote}}) as f, open({{out.outfile | quote}}, 'w') as fout:
	for line in f:
		line  = line.rstrip('\n')
		if not line: continue
		if line.startswith('#'):
			fout.write(line + '\n')
			continue
		parts  = line.split('\t')
		start  = int(parts[1])
		end    = int(parts[2])
		strand = parts[5]
		if not stdns or strand == '+':
			left2, right2 = right, left
		else:
			left2, right2 = left, right
		if 'pct' in params and params['pct']:
			length  = end - start
			start  -= round(length * left2)
			end    += round(length * right2)
		else:
			start -= left2
			end   += right2
		parts[1] = str(start)
		parts[2] = str(end)
		fout.write('\t'.join(parts) + '\n')
{% else %}
cmd = '{{args.bedtools}} flank %s > {{out.outfile | quote}}' % params2CmdArgs(params, dash='-', equal=' ')
runcmd(cmd)
{% endif %}