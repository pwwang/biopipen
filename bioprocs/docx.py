from pyppl import Proc, Box
from . import params

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
	`code`: Some extra code
@requires:
	[`python-docx`](http://python-docx.readthedocs.io/en/latest)
"""
pDocx             = Proc(desc = 'Operating .docx file')
pDocx.input       = 'infile, codes:files'
pDocx.output      = 'outfile:file:{{in.infile | lambda f, path = __import__("os").path: bn(f) if path.isfile(f) else str2fn(f) + ".docx"}}'
pDocx.args.code   = []
pDocx.envs.str2fn = lambda s, re = __import__('re'): re.sub(r'[^\w\-_\.]', '_', s)[:32]
pDocx.lang        = params.python.value
pDocx.script      = "file:scripts/docx/pDocx.py"

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
	`style`: The style of the table
	`bcode`: Some extra BEFORE the table is inserted.
	`acode`: Some extra AFTER the table is inserted.
@requires:
	[`python-docx`](http://python-docx.readthedocs.io/en/latest)
"""
pTable2DocxCode             = Proc(desc = 'Convert a table to docx code')
pTable2DocxCode.input       = 'infile:file'
pTable2DocxCode.output      = 'outfile:file:{{in.infile | fn2}}.docxcode.py'
pTable2DocxCode.args.style  = None
pTable2DocxCode.args.inopts = Box(rnames = True, cnames = True, skip = 0, delimit = "\t")
pTable2DocxCode.args.bcode  = []
pTable2DocxCode.args.acode  = []
pTable2DocxCode.lang        = params.Rscript.value
pTable2DocxCode.script      = "file:scripts/docx/pTable2DocxCode.r"

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
		- With keys `width` and/or `height` in inches.
	`bcode`: Some extra BEFORE the table is inserted.
	`acode`: Some extra AFTER the table is inserted.
@requires:
	[`python-docx`](http://python-docx.readthedocs.io/en/latest)
"""
pImage2DocxCode            = Proc(desc = 'Convert an image to docx code')
pImage2DocxCode.input      = 'infile:file'
pImage2DocxCode.output     = 'outfile:file:{{in.infile | fn2}}.docxcode.py'
pImage2DocxCode.args.scale = Box() # width and/or height in inches
pImage2DocxCode.args.bcode = []
pImage2DocxCode.args.acode = []
pImage2DocxCode.lang       = params.python.value
pImage2DocxCode.script     = "file:scripts/docx/pImage2DocxCode.py"