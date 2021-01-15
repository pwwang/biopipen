import re
import openpyxl
import csv

infile = {{ i.infile | quote}}
outdir = {{ o.outdir | quote}}
prefix = {{ o.outdir, stem(i.infile) | @join: '/' | quote}}

wb = openpyxl.load_workbook(infile, read_only = True)
for i, sh in enumerate(wb.sheetnames):
	sheet = wb.worksheets[i]
	outfile = prefix + '.' + re.sub(r'[^\w_]', '', sh) + '.tsv'

	with open(outfile, 'w', newline = "") as f:
		c = csv.writer(f, delimiter = '\t', escapechar = '\\', lineterminator = '\n',
			quoting = csv.QUOTE_NONE, doublequote = False)
		for r in sheet.rows:
			c.writerow([cell.value for cell in r])
