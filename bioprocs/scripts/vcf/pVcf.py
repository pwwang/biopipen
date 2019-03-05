import vcf
from pyppl import Box
from copy import copy

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
helper  = {{args.helper | repr}}
if not isinstance(helper, list):
	helper = [helper]
helper  = [line for line in helper if line]
exec('\n'.join(helper), globals())

reader    = vcf.Reader(filename=infile)
readerops = {{args.reader}}
newreader = readerops(copy(reader)) if callable(readerops) else None
writer    = vcf.Writer(open(outfile, 'w'), newreader if newreader else reader)
recordops = {{args.record}}

for record in reader:
	formats = record.FORMAT.split(':')
	samples = [
		Box(data = Box(zip(formats, [getattr(sample.data, fm) for fm in formats]))) 
		for sample in record.samples
	]
	rec     = recordops(record, writer, samples)
	for i, sample in enumerate(record.samples):
		sample.data = sample.data._make(list(samples[i].data.values()))
	if rec:
		writer.write_record(rec)
writer.close()
