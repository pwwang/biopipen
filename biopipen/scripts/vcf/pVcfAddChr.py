
from bioprocs.utils import FileConn
infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
chrom   = {{args.chr | quote}}

with FileConn(infile) as f, open(outfile, 'w') as fout:
	for line in f:
		if line.startswith('##contig='):
			contig = line.rstrip('>\n')[10:].split(',')
			items  = []
			for cont in contig:
				name, val = cont.split('=', 1)
				if name == 'ID':
					val = chrom + val if not val.startswith(chrom) else val
				items.append('{}={}'.format(name, val))
			fout.write('##contig=<{}>\n'.format(','.join(items)))
		elif line.startswith('#') or line.startswith(chrom):
			fout.write(line)
		else:
			fout.write(chrom + line)
