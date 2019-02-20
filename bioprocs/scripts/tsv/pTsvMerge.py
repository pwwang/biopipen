from pyppl import Box
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

files   = {{i.infiles | repr}}
inopts  = {{args.inopts | repr}}
outopts = {{args.outopts | repr}}
outfile = {{o.outfile | quote}}
maxopen = {{args.maxopen | repr}}

def getReaders(fs, baseidx):
	readers = []
	for i,infile in enumerate(fs):
		reader = TsvReader(infile, **inopts)
		if not reader.meta:
			reader.autoMeta()
		readers.append(reader)
	return readers

writer = TsvWriter(outfile, delimit = inopts['delimit'])
if inopts.get('cnames', True):
	reader = TsvReader(files[0], **inopts)
	reader.close()
	writer.cnames = reader.cnames
	writer.writeHead()

for i in range(0, len(files), maxopen):
	fs = files[i:i + maxopen]
	readers = getReaders(fs, i)

	for reader in readers:
		for r in reader:
			writer.write(r)
		reader.close()
writer.close()
