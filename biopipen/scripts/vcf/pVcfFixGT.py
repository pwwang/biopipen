from bioprocs.utils.shell2 import bgzip

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
gz      = {{args.gz | bool}}
if gz:
	outfile = outfile[:-3]

with open(infile) as fin, open(outfile, 'w') as fout:
	for line in fin:
		line = line.rstrip('\n')
		if line.startswith('#'):
			fout.write(line + '\n')
			continue

		parts = line.split('\t')
		for i, part in enumerate(parts):
			formats = part.split(':')
			if formats[0] != '0/0' and formats[0] != '0|0':
				continue
			if any(fmt != '.' for fmt in formats[1:]):
				continue
			formats[0] = './.'
			parts[i] = ':'.join(formats)
		fout.write('\t'.join(parts) + '\n')

if gz:
	bgzip(outfile)
