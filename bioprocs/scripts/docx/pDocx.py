import docx
from docx.shared import Inches
from os import path
from runpy import run_path

infile    = {{in.infile | quote}}
outfile   = {{out.outfile | quote}}
codefiles = {{in.codes | repr}}

if path.isfile(infile):
	doc = docx.Document(infile)
else:
	doc = docx.Document()
	doc.add_heading(infile, 0)

for codefile in codefiles:
	_ = run_path(codefile, globals())

doc.save(outfile)



