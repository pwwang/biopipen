import openpyxl
import csv

infile = {{ i.infile | quote}}
outfile = {{ o.outfile | quote}}

wb = openpyxl.load_workbook(infile, read_only = True)
sh = wb.get_active_sheet()
with open(outfile, 'w', newline = "") as f:
    c = csv.writer(f, delimiter = '\t', escapechar = '\\', lineterminator = '\n',
		quoting = csv.QUOTE_NONE, doublequote = False)
    for r in sh.rows:
        c.writerow([cell.value for cell in r])
