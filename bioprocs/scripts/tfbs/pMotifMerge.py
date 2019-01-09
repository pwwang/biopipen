from bioprocs.utils.meme import MemeReader, MemeWriter

infiles = {{ i.infiles | repr}}
outfile = {{ o.outfile | quote}}

reader0 = MemeReader(infiles[0])
writer  = MemeWriter(outfile)
writer.meta = reader0.meta
writer.writeMeta()
reader0.close()

for infile in infiles:
	reader = MemeReader(infile)
	for r in reader:
		writer.write(r)
	reader.close()

writer.close()

