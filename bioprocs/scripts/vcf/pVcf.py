import vcf
from copy import copy

infile    = {{in.infile | quote}}
{{args.helper}}


reader    = vcf.Reader(filename=infile)
readerops = {{args.readerops}}
newreader = readerops(copy(reader))
writer    = vcf.Writer(open({{out.outfile | quote}}, 'w'), newreader if newreader else reader)
recordops = {{args.recordops}}
for record in reader:
	rec = recordops(record, writer)
	if rec:
		writer.write_record(rec)
writer.close()
