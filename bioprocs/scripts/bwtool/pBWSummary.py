from pyppl import Box
from os import path
from bioprocs.utils import shell2 as shell
from bioprocs.utils.tsvio2 import TsvWriter
from gff import Gff

infile  = {{i.infile | quote}}
bedfile = {{i.bedfile | quote}}
bedext  = {{i.bedfile | ?:ext(_) == '.gz' | prefix | :_ | ext | [1:] | quote }}
outfile = {{o.outfile | quote}}
bwtool  = {{args.bwtool | quote}}
params  = {{args.params | repr}}

if 'bed' not in bedext:
	def getname(record):
		for key, val in record['attributes'].items():
			if key.endswith('id'):
				return val
		for key, val in record['attributes'].items():
			if key.endswith('name'):
				return val
		return f'{record["seqid"]}:{record["start"]}-{record["end"]}'

	tmpfile = path.join(path.dirname(outfile), 'regions.bed')
	gff = Gff(bedfile)
	writer = TsvWriter(tmpfile)
	for r in gff:
		writer.write([r['seqid'], r['start'], r['end'], getname(r)])
	writer.close()
	bedfile = tmpfile

shell.load_config(bwtool = dict(_exe = bwtool, _prefix = '-'))

shell.fg.bwtool.summary(bedfile, infile, outfile, **params)
