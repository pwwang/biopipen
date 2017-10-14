from pyppl import Proc, Box
from .utils import helpers, runcmd
from . import params

pTabix                        = Proc(desc = 'Use tabix to extract information.')
pTabix.input                  = "infile, region"
pTabix.output                 = "outfile:file:{{in.infile | fn | fn}}-{{job.index}}{% if in.infile.endswith('.gz') %}{{in.infile | fn | ext}}{% else %}{{in.infile | ext}}{% endif %}"
pTabix.args.tabix             = params.tabix.value
pTabix.args.params            = Box({'h': True})
pTabix.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pTabix.tplenvs.runcmd         = runcmd.py
pTabix.lang                   = params.python.value
pTabix.script                 = "file:scripts/tabix/pTabix.py"

pTabixIndex                        = Proc(desc = 'Generate tabix index file')
pTabixIndex.input                  = "infile:file"
pTabixIndex.output                 = [
	"outfile:file:{{in.infile | bn}}{% if in.infile.endswith('.gz') | lambda x: not x %}.gz{% endif %}", 
	"outidx:file:{{in.infile | bn}}{% if in.infile.endswith('.gz') | lambda x: not x %}.gz{% endif %}.tbi"
]
pTabixIndex.args.tabix             = params.tabix.value
pTabixIndex.args.python            = params.python.value
pTabixIndex.args.params            = Box()
pTabixIndex.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pTabixIndex.script                 = "file:scripts/tabix/pTabixIndex.bash"