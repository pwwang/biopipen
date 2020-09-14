
import csv

infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}

with open(infile, newline = '') as csvfile:
	sample = csvfile.read(10240)
	csvfile.seek(0)
	dialect = csv.Sniffer().sniff(sample)
	has_header = csv.Sniffer().has_header(sample)
	reader = csv.reader(csvfile, dialect)
	with open(outfile, 'w', newline = "") as f:
		cout = csv.writer(f, delimiter = '\t', escapechar = '\\', lineterminator = '\n',
			quoting = csv.QUOTE_NONE, doublequote = False)
		cout.writerows(reader)
