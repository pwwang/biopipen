import docx
from docx.shared import Inches, Pt
from docx.enum.dml import MSO_COLOR_TYPE, MSO_THEME_COLOR
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING, WD_BREAK_TYPE, WD_COLOR, WD_TAB_ALIGNMENT, WD_TAB_LEADER, WD_UNDERLINE
from docx.enum.style import WD_STYLE, WD_STYLE_TYPE
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_TABLE_DIRECTION
from docx.enum.shape import WD_INLINE_SHAPE_TYPE
from os import path
from runpy import run_path
from bioprocs.utils import log2pyppl

def secOri(s, ori):
	s.orientation = ori
	new_width, new_height = s.page_height, s.page_width
	s.page_width  = new_width
	s.page_height = new_height

infile    = {{in.infile | quote}}
outfile   = {{out.outfile | quote}}
codefiles = {{in.codes | repr}}
error     = {{args.error | repr}}
bcode     = {{args.bcode | repr}}
if not isinstance(bcode, list):
	bcode = [bcode]
bcode = '\n'.join(bcode) + '\n'
acode     = {{args.acode | repr}}
if not isinstance(acode, list):
	acode = [acode]
acode = '\n'.join(acode) + '\n'

if path.isfile(infile):
	doc = docx.Document(infile)
	exec(bcode, globals())
else:
	doc = docx.Document()
	exec(bcode, globals())
	doc.add_heading(infile, 0)


for codefile in codefiles:
	try:
		_ = run_path(codefile, globals())
	except Exception as ex:
		if error == 'exit':
			raise
		else:
			log2pyppl('Failed to run: {}'.format(codefile))
			for line in str(ex).splitlines():
				log2pyppl('\t' + line)

exec(acode, globals())
doc.save(outfile)



