"""A set of processes handling Excel files (>=2010)"""

from pyppl import Proc
from . import params
from .utils import fs2name
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pTsvs2Xlsx():
	"""
	@name:
		pTsvs2Xlsx
	@description:
		Save tsv files to xlsx sheets.
	@input:
		`infiles:files`: The input tsv files
	@output:
		`outfile:file`: The output xlsx file
	@args:
		`fn2sheet`: How to convert filename(without extension) to sheet name
	@requires:
		python packages: `csv` and `openpyxl`
	"""
	pTsvs2Xlsx               = Proc(desc = 'Save tsv files to xlsx sheets.')
	pTsvs2Xlsx.input         = 'infiles:files'
	pTsvs2Xlsx.output        = 'outfile:file:{{i.infiles | fs2name}}.xlsx'
	pTsvs2Xlsx.args.fn2sheet = None
	pTsvs2Xlsx.envs.fs2name  = fs2name
	pTsvs2Xlsx.lang          = params.python.value
	pTsvs2Xlsx.script        = "file:scripts/xlsx/pTsvs2Xlsx.py"
	return pTsvs2Xlsx

