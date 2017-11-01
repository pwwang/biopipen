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
		if 'pct' in params and params['pct']:
			length  = end - start
			start  -= round(length * left)
			end    += round(length * right)
		else:
			start -= left
			end   += right
		parts[1] = str(start)
		parts[2] = str(end)
		fout.write('\t'.join(parts) + '\n')
{% else %}
cmd = '{{args.bedtools}} flank %s > {{out.outfile | quote}}' % params2CmdArgs(params, dash='-', equal=' ')
runcmd(cmd)
{% endif %}