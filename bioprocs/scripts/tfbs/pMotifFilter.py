from bioprocs.utils.meme import MemeReader, MemeWriter

infile  = {{ i.infile | quote}}
outfile = {{ o.outfile | quote}}
# if filter has multiple lines, treat the first lines as helper and last line as the funciton
# now load the first lines
{{ args.filter | :"" if not "\n" in a else "\n".join(a.split("\n")[:-1]) }}
ffunc   = {{ args.filter | :None if a is None else a.split("\n")[-1] }}

reader = MemeReader(infile)
writer = MemeWriter(outfile)
writer.meta = reader.meta
writer.writeMeta()

for r in reader:
	if callable(ffunc) and ffunc(r):
		continue
	writer.write(r)
reader.close()
writer.close()
