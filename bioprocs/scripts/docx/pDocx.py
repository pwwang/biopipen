import docx
from pyppl import Box
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

infile    = {{i.infile | quote}}
outfile   = {{o.outfile | quote}}
codefiles = {{i.codes | repr}}
error     = {{args.error | repr}}
bcode     = {{args.bcode | repr}}
if not isinstance(bcode, list):
	bcode = [bcode]
bcode = '\n'.join(bcode) + '\n'
acode     = {{args.acode | repr}}
if not isinstance(acode, list):
	acode = [acode]
acode = '\n'.join(acode) + '\n'
section   = {{args.section | repr}}
if section:
	section['type']   = section.get('type', 'NEW_PAGE')
	section['orient'] = section.get('orient', 'PORTRAIT')
	section['margin'] = section.get('margin', [36] * 4)
	if isinstance(section['margin'], int):
		section['margin'] = [section['margin']] * 4
	elif len(section['margin']) == 2:
		section['margin'].extend(section['margin'])
	elif len(section['margin']) == 3:
		section['margin'].append(section['margin'][1])

def doSection(doc, section = section, new = True):
	if not section:
		return
	sec = doc.add_section() if new else doc.sections[0]
	sec.type = getattr(WD_SECTION, section['type'])
	sec.orientation = getattr(WD_ORIENT, section['orient'])
	_new_width, _new_height = sec.page_height, sec.page_width
	sec.page_width, sec.page_height = _new_width, _new_height
	sec.top_margin, sec.right_margin, sec.bottom_margin, sec.left_margin = (Pt(x) for x in section['margin'])
	return sec

if infile and path.isfile(infile):
	doc = docx.Document(infile)
	# do section
	sec = doSection(doc, new = True)
	exec(bcode, globals())
else:
	doc = docx.Document()
	sec = doSection(doc, new = False)
	if infile:
		doc.add_heading(infile, 0)
	exec(bcode, globals())

for codefile in codefiles:
	log2pyppl('Doing: {}'.format(path.basename(codefile)), level = 'Info')
	try:
		_ = run_path(codefile, globals())
	except Exception as ex:
		if error == 'exit':
			raise
		else:
			log2pyppl('Failed to run: {}'.format(codefile), level = 'Error')
			for line in str(ex).splitlines():
				log2pyppl('\t' + line)

exec(acode, globals())
doc.save(outfile)



