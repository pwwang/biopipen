import random
import vcf

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
keep    = {{args.keep | quote}}
upart   = {{args.upart | repr}}
gz      = {{args.gz | repr}}
if gz:
	outfile = outfile[:-3]

keeps = [k.strip() for k in keep.split(',')]

reader = vcf.Reader(filename=infile)
writer = vcf.Writer(open(outfile, 'w'), reader)

def writeRecord(records, ks):
	if len(records) == 1:
		writer.write_record(records[0])
	else:
		if 'bisnp' in ks:
			writeRecord([r for r in records if r.is_snp and len(r.ALT) == 1], [k for k in ks if k!='bisnp'])
		elif 'snp' in ks:
			writeRecord([r for r in records if r.is_snp], [k for k in ks if k!='snp'])
		elif 'bialt' in ks:
			writeRecord([r for r in records if len(r.ALT) == 1], [k for k in ks if k!='bialt'])
		elif 'last' in ks:
			writer.write_record(records[-1])
		elif 'random' in ks:
			writer.write_record(random.choice(records))
		else:
			writer.write_record(records[0])

buffer_name = None
buffer_recs = []
while True:
	try:
		record = reader.next()
		rid = '/'.join([str(getattr(record, i)) for i in upart])
		if buffer_name is None:
			buffer_name = rid
			buffer_recs.append(record)
		elif buffer_name == rid:
			buffer_recs.append(record)
		else:
			writeRecord(buffer_recs, keeps)
			buffer_name = rid
			buffer_recs = [record]
	except StopIteration:
		break
	except ValueError:
		continue

writeRecord(buffer_recs, keeps)
writer.close()

if gz:
	runcmd(['bgzip', outfile])


