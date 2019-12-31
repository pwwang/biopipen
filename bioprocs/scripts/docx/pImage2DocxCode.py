from diot import Diot

bcode   = {{args.bcode | repr}}
acode   = {{args.acode | repr}}
infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
scale   = {{args.scale | repr}}
section = {{args.section | repr}}
legend  = {{args.legend | render | repr}}
{% case True %}
{% when isinstance(args.title, (tuple, list)) %}
title = {{args.title | lambda a, render = render: (render(a[0]), a[1])}}
{% when args.title | bool %}
title = {{args.title | lambda a, render = render: (render(a), 1)}}
{% else %}
title = {{args.title | repr}}
{% endcase %}
align   = {{args.align | repr}}

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

scale = ', '.join('{} = Pt({})'.format(k, v) for k, v in scale.items())
scale = ', ' + scale if scale else ''

if legend:
	legend = legend.split('.', 1) if legend.startswith('Figure ') and '.' in legend else ['', legend]

# region template
template = """
# Set up section
{section_code}

# Add heading
{heading_code}

# Interpolation code before image
{bcode_code}

# Insert the image
{image_code}

# Add legend
{legend_code}

# Interpolation code after image
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

image_code = """
para = doc.add_paragraph()
{alignment}
run = para.add_run()
run.add_picture({pic!r}{scale})
""".format(
	alignment = "para.alignment = WD_ALIGN_PARAGRAPH.{}\n".format(align.upper()) if align else '',
	pic       = infile,
	scale     = scale
)

legend_code = """
_legend_para = doc.add_paragraph()
{alignment}
_legend_para.add_run({legend0!r}).font.bold = True
_legend_para.add_run({legend1!r})
""".format(
	alignment = "_legend_para.alignment = WD_ALIGN_PARAGRAPH.{}\n".format(align.upper()) if align else '',
	legend0   = legend[0],
	legend1   = legend[1]
) if legend else ''

acode_code = acode

with open(outfile, 'w') as f:
	f.write(template.format(
		section_code = section_code,
		heading_code = heading_code,
		bcode_code   = bcode_code,
		image_code   = image_code,
		legend_code  = legend_code,
		acode_code   = acode_code
	))
