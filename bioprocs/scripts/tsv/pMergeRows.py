from pyppl import Box
from bioprocs.utils.tsvio import TsvReader, TsvWriter, TsvRecord

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
inopts  = {{args.inopts}}
outopts = {{args.outopts}}
match   = {{args.match}}
do      = {{args.do}}

reader = TsvReader(infile, **inopts)
writer = TsvWriter(outfile, **outopts)

if not writer.meta:
	writer.meta = reader.meta

if outopts.head:
	writer.writeHead(
		delimit   = outopts.headDelimit,
		prefix    = outopts.headPrefix,
		transform = outopts.headTransform
	)

prev     = None
repeated = []
for r in reader:
	m = r[0] if not callable(match) else match(r)
	if prev is None:
		prev = m
		repeated.append(r)
	elif prev == m:
		repeated.append(r)
	else:
		prev = m
		writer.write(do(repeated))
		repeated = [r]

if repeated:
	writer.write(do(repeated))

	