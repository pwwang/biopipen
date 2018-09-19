from os import path, makedirs, remove
from shutil import rmtree, copyfile, move
from sys import stderr
import vcf
from vcf.model import _Call
from bioprocs.utils import runcmd

filters = {{args.filters}}
gz      = {{args.gz | repr}}
outfile = {{o.outfile | quote}}
if gz:
	outfile = outfile[:-3]

reader = vcf.Reader(filename='{{i.infile}}')
writer = vcf.Writer(open(outfile, 'w'), reader)
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

if gz:
	runcmd(['bgzip', outfile])


