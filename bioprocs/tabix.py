from pyppl import Proc, Box
#from .utils import helpers, runcmd
from . import params

"""
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
"""
pTabix             = Proc(desc = 'Use tabix to extract information.')
pTabix.input       = "infile, region"
pTabix.output      = "outfile:file:{{in.infile | fn | fn}}-{{job.index}}{% if in.infile.endswith('.gz') %}{{in.infile | fn | ext}}{% else %}{{in.infile | ext}}{% endif %}"
pTabix.args.tabix  = params.tabix.value
pTabix.args.params = Box(h = True)
pTabix.lang        = params.python.value
pTabix.script      = "file:scripts/tabix/pTabix.py"

"""
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
	`python`: Will be used to generate command line arguments.
"""
pTabixIndex        = Proc(desc = 'Generate tabix index file')
pTabixIndex.input  = "infile:file"
pTabixIndex.output = [
	"outfile:file:{{in.infile | bn}}{% if in.infile.endswith('.gz') | lambda x: not x %}.gz{% endif %}", 
	"outidx:file:{{in.infile | bn}}{% if in.infile.endswith('.gz') | lambda x: not x %}.gz{% endif %}.tbi"
]
pTabixIndex.args.tabix  = params.tabix.value
pTabixIndex.args.python = params.python.value
pTabixIndex.args.params = Box()
pTabixIndex.script      = "file:scripts/tabix/pTabixIndex.bash"