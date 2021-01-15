import csv
from os import path
from openpyxl import Workbook

infile  = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
fn2sheet = {{args.fn2sheet}}

def tsv2sheet(wb, tsvfile):
	ws = wb.create_sheet('Sheet1')
	with open(tsvfile) as f:
		reader = csv.reader(f, delimiter = "\t")
		for row in reader:
			ws.append(row)

wb = Workbook()
wb.remove(wb.active) # remove default sheet
tsv2sheet(wb, infile)

wb.save(outfile)
