from pyppl import Box
from bioprocs.utils.tsvio2 import TsvReader

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
font    = {{args.font | repr}}
section = {{args.section | repr}}
{% case True %}
{% when isinstance(args.title, (tuple, list)) %}
title = {{args.title | lambda a, render = render: (render(a[0]), a[1])}}
{% when args.title | bool %}
title = {{args.title | lambda a, render = render: (render(a), 1)}}
{% else %}
title = {{args.title | repr}}
{% endcase %}
caption = {{args.caption | render | repr}}
style   = {{args.style | repr}}
inopts  = {{args.inopts | repr}}
acode   = {{args.acode | repr}}
bcode   = {{args.bcode | repr}}

# fix parameters
if section:
	section['type']   = section.get('type', 'NEW_PAGE').upper()
	section['orient'] = section.get('orient', 'PORTRAIT').upper()
	section['margin'] = section.get('margin', [36] * 4)
	if isinstance(section['margin'], int):
		section['margin'] = [section['margin']] * 4
	elif len(section['margin']) == 2:
		section['margin'].extend(section['margin'])
	elif len(section['margin']) == 3:
		section['margin'].append(section['margin'][1])

if not isinstance(bcode, (tuple, list)):
	bcode = [bcode]
bcode = "\n".join(bcode)

if not isinstance(acode, (tuple, list)):
	acode = [acode]
acode = "\n".join(acode)

if caption:
	caption = caption.split('.', 1) if caption.startswith('Table ') and '.' in caption else ['', caption]

# region template
template = """
# Set up section
{section_code}

# Add heading
{heading_code}

# Interpolation code before table
{bcode_code}

# Add caption
{caption_code}

# Insert the table
{table_code}

# Interpolation code after table
{acode_code}
"""
# endregion

# region section_code
section_code = """
sec             = doc.add_section()
sec.type        = WD_SECTION.{sectype}
sec.orientation = WD_ORIENT.{secorient}

_new_width, _new_height = sec.page_height, sec.page_width
if sec.orientation != doc.sections[-2].orientation:
	sec.page_width, sec.page_height = _new_width, _new_height
else:
	sec.page_width, sec.page_height = _new_height, _new_width
sec.top_margin, sec.right_margin, sec.bottom_margin, sec.left_margin = (Pt(x) for x in {secmargin!r})
""".format(
	sectype   = section['type'],
	secorient = section['orient'],
	secmargin = section['margin']
) if section else ""
# endregion

# region heading_code
heading_code = """
doc.add_heading('''{}''', {})
""".format(*title) if title else ""
# endregion

bcode_code = bcode

caption_code = """
_cap_para = doc.add_paragraph()
{alignment}
_cap_para.add_run({caption0!r}).font.bold = True
_cap_para.add_run({caption1!r})
""".format(
	alignment = "_cap_para.alignment = WD_ALIGN_PARAGRAPH.{}\n".format(align.upper()) if align else '',
	caption0   = caption[0],
	caption1   = caption[1]
) if caption else ''

inopts['cname0'] = inopts.get('cname0', '')
reader = TsvReader(infile, **inopts)
table_code = """
cnames = {cnames!r}
data   = {data!r}
if cnames:
	data.insert(0, cnames)
table = doc.add_table(rows = len(data), cols = len(cnames))
table.style = {style!r}
for i, row in enumerate(data):
	for j, r in enumerate(row):
		table.rows[i].cells[j].text = r
		table.rows[i].cells[j].paragraphs[0].runs[0].font.name = {family!r}
		table.rows[i].cells[j].paragraphs[0].runs[0].font.size = Pt({size!r})
""".format(
	cnames = reader.cnames,
	data   = [r.values() for r in reader],
	family = font['family'],
	size   = font['size'],
	style  = style
)

acode_code = acode

with open(outfile, 'w') as f:
	f.write(template.format(
		section_code = section_code,
		heading_code = heading_code,
		bcode_code   = bcode_code,
		caption_code = caption_code,
		table_code   = table_code,
		acode_code   = acode_code
	))

