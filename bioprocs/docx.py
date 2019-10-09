"""Processes for WORD files"""
from pyppl import Proc, Box
from . import params, rimport
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pDocx():
	"""
	@name:
		pDocx
	@description:
		Operating .docx file
	@input:
		`infile`: A input docx file or a string that used as Heading 0 of a new file
		`codes:files`: A set of files containing python codes operating the docs file.
			- Variable `doc` is used to pass across files.
	@output:
		`outfile:file`: The output docx file
	@args:
		`bcode`: Some extra BEFORE the content is inserted (after heading).
		`acode`: Some extra AFTER the content is inserted.
		`error`: What to do when error happens. Default: `exit`
			- `ignore` to add nothing to the document
		`section`: Start a new section of this table. Could be a box of values. Default: `Box()`.
			- `Box()` means no new section will start.
			- `type`: The type of the section. Could be one of:
				- `CONTINUOUS`
				- `NEW_COLUMN`
				- `NEW_PAGE` (default if not provided when `args.section.type` is not `None`)
				- `EVEN_PAGE` and
				- `ODD_PAGE`
			- `orient`: The orientation of the section. Default: `PORTRAIT`
				- Could also be `LANDSCAPE`
			- `margin`: The margin of the section, in points. Could be:
				- Single value for all margins. Default: 36
				- Paired values for top/bottom and left/right margins
				- 3 values for top, left/right and bottom margins.
				- 4 values for top, right, bottom and left margins
	@requires:
		[`python-docx`](http://python-docx.readthedocs.io/en/latest)
	"""
	pDocx              = Proc(desc = 'Operating .docx file')
	pDocx.input        = 'infile, codes:files'
	pDocx.output       = 'outfile:file:{{bn(i.infile) if path.isfile(i.infile) else str2fn(i.infile) + ".docx"}}'
	pDocx.args.bcode   = []
	pDocx.args.acode   = []
	pDocx.args.error   = 'exit' # or ignore
	pDocx.args.section = Box()
	pDocx.envs.str2fn  = lambda s, re = __import__('re'): re.sub(r'[^\w\-_\.]', '_', s)[:32] if s else 'doc'
	pDocx.envs.path    = __import__("os").path
	pDocx.lang         = params.python.value
	pDocx.script       = "file:scripts/docx/pDocx.py"
	return pDocx

@procfactory
def _pTable2DocxCode():
	"""
	@name:
		pTable2DocxCode
	@description:
		Convert a table to docx code
	@input:
		`infile`: The table data file
	@output:
		`outfile:file`: The code file
	@args:
		`font`   : The font of table content. Default: `Box(family: 'Arial', size: 10)`
		`section`: Start a new section of this table. Could be a box of values. Default: `Box()`.
			- `Box()` means no new section will start.
			- `type`: The type of the section. Could be one of:
				- `CONTINUOUS`
				- `NEW_COLUMN`
				- `NEW_PAGE` (default if not provided when `args.section.type` is not `None`)
				- `EVEN_PAGE` and
				- `ODD_PAGE`
			- `orient`: The orientation of the section. Default: `PORTRAIT`
				- Could also be `LANDSCAPE`
			- `margin`: The margin of the section, in points. Could be:
				- Single value for all margins. Default: 36
				- Paired values for top/bottom and left/right margins
				- 3 values for top, left/right and bottom margins.
				- 4 values for top, right, bottom and left margins
		`style`  : The style of the table. Default: `Light List Accent 1`
			- `Table Normal`
			- `Colorful Grid`, `Colorful Grid Accent 1~6`
			- `Colorful List`, `Colorful List Accent 1~6`
			- `Colorful Shading`, `Colorful Shading Accent 1~6`
			- `Dark List`, `Dark List Accent 1~6`
			- `Light Grid`, `Light Grid Accent 1~6`
			- `Light List`, `Light List Accent 1~6`
			- `Light Shading`, `Light Shading Accent 1~6`
			- `Medium Grid 1~3`, `Medium Grid 1~3 Accent 1~6`
			- `Medium List 1~2`, `Medium List 1~2 Accent 1~6`
			- `Medium Shading 1~2`, `Medium Shading 1~2 Accent 1~6`
			- `Table Grid`
		`title`  : The title of the table, could be heading to table caption, Could be:
			- A string, defaulted as `Heading 1` (2nd level heading)
			- A `tuple` or `list` with first element (string) as the content, 2nd element as the heading level
		`caption`: The caption of the table. Default: `None`
		`bcode`  : Some extra BEFORE the table is inserted.
		`acode`  : Some extra AFTER the table is inserted.
		`inopts` : The options for reading the input file.
	@requires:
		[`python-docx` v0.8.7](http://python-docx.readthedocs.io/en/latest)
	"""
	pTable2DocxCode              = Proc(desc = 'Convert a table to docx code')
	pTable2DocxCode.input        = 'infile:file'
	pTable2DocxCode.output       = 'outfile:file:{{i.infile | fn2}}.docxcode.py'
	pTable2DocxCode.args.font    = Box(family = 'Arial', size = 10)
	pTable2DocxCode.args.section = Box()
	pTable2DocxCode.args.style   = 'Light List Accent 1'
	pTable2DocxCode.args.title   = None
	pTable2DocxCode.args.caption = None
	pTable2DocxCode.args.inopts  = Box(cnames = True, skip = 0, delimit = "\t")
	pTable2DocxCode.args.bcode   = []
	pTable2DocxCode.args.acode   = []
	pTable2DocxCode.envs.rimport = rimport
	pTable2DocxCode.lang         = params.python.value
	pTable2DocxCode.script       = "file:scripts/docx/pTable2DocxCode.py"
	return pTable2DocxCode

@procfactory
def _pImage2DocxCode():
	"""
	@name:
		pImage2DocxCode
	@description:
		Convert an image to docx code
	@input:
		`infile`: The image
	@output:
		`outfile:file`: The code file
	@args:
		`scale`: The scale of the image. Default: `Box()`
			- With keys `width` and/or `height` in points.
			- If only one is set, the other is automatically scaled.
		`title`  : The title of the table, could be heading to table caption, Could be:
			- A string, defaulted as `Heading 1` (2nd level heading)
			- A `tuple` or `list` with first element (string) as the content, 2nd element as the heading level
			- Template available. i.e. `("Title: {{i.infile | fn2}}", 2)`
		`section`: Start a new section of this table. Could be a box of values. Default: `Box()`.
			- `Box()` means no new section will start.
			- `type`: The type of the section. Could be one of:
				- `CONTINUOUS`
				- `NEW_COLUMN`
				- `NEW_PAGE` (default if not provided when `args.section.type` is not `None`)
				- `EVEN_PAGE` and
				- `ODD_PAGE`
			- `orient`: The orientation of the section. Default: `PORTRAIT`
				- Could also be `LANDSCAPE`
			- `margin`: The margin of the section, in points. Could be:
				- Single value for all margins. Default: 36
				- Paired values for top/bottom and left/right margins
				- 3 values for top, left/right and bottom margins.
				- 4 values for top, right, bottom and left margins
		`legend`: The legend of the image. Default: `None`
			- If starts with `Figure X. `, this part will be bolded.
			- Template available. i.e. `Figure 1. {{i.infile | fn2}} shows good results.`
		`align `: How the image is aligned if it is not an inline element. Default: `CENTER`
			- To set it inline, set `args.align` as `None`, and add something to the paragraph(`para`) in `args.acode`.
		`bcode` : Some extra BEFORE the table is inserted.
			- `doc` is available for the `document` object
			- `sec` is available for the `section` object if `args.section` has been set.
		`acode` : Some extra AFTER the table is inserted.
			- `doc` and `sec` are still available
			- `para`: The paragraph where the run of the image is placed
			- `run`: The run where the image is placed
	@requires:
		[`python-docx`](http://python-docx.readthedocs.io/en/latest)
	"""
	pImage2DocxCode              = Proc(desc = 'Convert an image to docx code')
	pImage2DocxCode.input        = 'infile:file'
	pImage2DocxCode.output       = 'outfile:file:{{i.infile | fn2}}.docxcode.py'
	pImage2DocxCode.args.scale   = Box() # width and/or height in inches
	pImage2DocxCode.args.section = Box()
	pImage2DocxCode.args.title   = None
	pImage2DocxCode.args.legend  = None
	pImage2DocxCode.args.align   = 'CENTER'
	pImage2DocxCode.args.bcode   = []
	pImage2DocxCode.args.acode   = []
	pImage2DocxCode.lang         = params.python.value
	pImage2DocxCode.script       = "file:scripts/docx/pImage2DocxCode.py"
	return pImage2DocxCode

