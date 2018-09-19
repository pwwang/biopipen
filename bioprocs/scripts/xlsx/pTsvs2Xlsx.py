import csv
from os import path
from openpyxl import Workbook
from bioprocs.utils import log2pyppl

infiles  = {{i.infiles | repr}}
outfile  = {{o.outfile | quote}}
fn2sheet = {{args.fn2sheet}}

def tsv2sheet(wb, tsvfile):
	fn = path.splitext(path.basename(tsvfile))[0]
	ws = wb.create_sheet(fn if not callable(fn2sheet) else fn2sheet(fn))
	log2pyppl('Reading ' + tsvfile + ' ...')
	with open(tsvfile) as f:
		reader = csv.reader(f, delimiter = "\t")
		for row in reader:
			ws.append(row)

wb = Workbook()
wb.remove(wb.active) # remove default sheet
for infile in infiles:
	tsv2sheet(wb, infile)

wb.save(outfile)