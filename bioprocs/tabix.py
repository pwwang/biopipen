"""Tabix utilities"""
from pyppl import Proc, Diot
#from .utils import helpers, runcmd
from . import params, proc_factory

pTabix = proc_factory(
	desc = 'Use tabix to extract information.',
	config = Diot(annotate = """
	@name:
		pTabix
	@description:
		Use tabix to extract information.
	@input:
		`infile`: a local or remote file
		`region`: a region or a file containing regions
	@output:
		`outfile:file`: The information extracted from the input file
	@args:
		`tabix`: The path to `tabix`
		`params`: Other params for `tabix`
	"""))
pTabix.input       = "infile, region"
pTabix.output      = "outfile:file:{{i.infile | fn | fn}}-{{job.index}}{% if i.infile.endswith('.gz') %}{{i.infile | fn | ext}}{% else %}{{i.infile | ext}}{% endif %}"
pTabix.args.tabix  = params.tabix.value
pTabix.args.params = Diot(h = True)
pTabix.lang        = params.python.value

pTabixIndex = proc_factory(
	desc = 'Generate tabix index file',
	config = Diot(annotate = """
	@name:
		pTabixIndex
	@description:
		Generate tabix index file.
	@input:
		`infile:file`: the input file
			- Could be bgzipped.
	@output:
		`outfile:file`: The bgzipped file
		`outidx:file`:  The tabix index file
	@args:
		`tabix`: The path to `tabix`
		`params`: Other params for `tabix`
	"""))
pTabixIndex.input  = "infile:file"
pTabixIndex.output = [
	"outfile:file:{{i.infile | bn}}{% if args.gz %}.gz{% endif %}",
	"outidx:file:{{i.infile  | bn}}{% if args.gz %}.gz{% endif %}.tbi"
]
pTabixIndex.args.gz     = True
pTabixIndex.args.tabix  = params.tabix.value
pTabixIndex.args.params = Diot()
pTabixIndex.lang        = params.python.value
