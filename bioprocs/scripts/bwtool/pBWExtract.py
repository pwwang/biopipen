from bioprocs.utils import shell2 as shell
from bioprocs.utils.tsvio2 import TsvWriter
from gff import Gff

infile  = {{i.infile | quote}}
bedfile = {{i.bedfile | quote}}
bedext  = {{i.bedfile | ?:ext(_) == '.gz' | prefix | :_ | ext | [1:] | quote }}
outfile = {{o.outfile | quote}}
bwtool  = {{args.bwtool | quote}}

if 'bed' not in bedext:
	tmpfile = path.join(path.dirname(outfile), 'regions.bed')
	gff = Gff(bedfile)
	writer = TsvWriter(tmpfile)
	for r in gff:
		writer.write([
			r['seqid'], r['start'], r['end'],
			r['attributes'].get(
				'id',
				r['attributes'].get('name', f'{r["seqid"]}:{r["start"]}-{r["end"]}'))
		])
	bedfile = tmpfile

shell.load_config(bwtool = dict(_exe = bwtool, _prefix = '-'))

shell.fg.bwtool.extract('bed', bedfile, infile, outfile, tabs = True)
