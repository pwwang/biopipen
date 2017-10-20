from os import path, makedirs, remove
from shutil import rmtree, copyfile, move
from sys import stderr
import vcf
from vcf.model import _Call


filters   = {{args.filters}}

reader = vcf.Reader(filename='{{in.infile}}')
writer = vcf.Writer(open('{{out.outfile}}', 'w'), reader)
for record in reader:
	if filters(record, record.samples):
		record.FILTER = ['PASS']
		writer.write_record(record)
	else:
		{% if args.keep %}
		record.FILTER  = ['{{args.fname}}']
		writer.write_record(record)
		{% else %}
		pass
		{% endif %}
writer.close()
