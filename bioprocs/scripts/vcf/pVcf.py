import vcf
from copy import copy

infile    = {{i.infile | quote}}
{{args.helper}}


reader    = vcf.Reader(filename=infile)
readerops = {{args.readerops}}
newreader = readerops(copy(reader))
writer    = vcf.Writer(open({{o.outfile | quote}}, 'w'), newreader if newreader else reader)
recordops = {{args.recordops}}
for record in reader:
	rec = recordops(record, writer)
	if rec:
		writer.write_record(rec)
writer.close()
