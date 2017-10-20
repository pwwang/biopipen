import re
from pyppl import Proc, Box
from .utils import runcmd, helpers, parallel
from . import params
"""
@name:
	pMotifScan
@input:
	`mfile:file`: The motif file
	`sfile:file`: The sequence file
@output:
	`outdir:file`: The output dir
@args:
	`tools`   : The tool used to scan the motif. Default: 'meme'
	`meme`    : The path of MEME's fimo. Default: 'fimo'
	`mids`    : The motif ids from the motif input file, default: "all"
		- Can be a list of multiple ids or a string separated by comma
	`pval`    : The pvalue cutoff. Default: 1e-4
	`cleanmname`: Whether to clean motif name. Default: True
	`ucsclink`: The ucsc link template. Default: `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={}`
	`nthread` : Number of threads used to scan, only available when you have multiple mids. Default: 1
	`params`  : Other parameters for `fimo`
@requires:
	[`fimo` from MEME Suite](http://meme-suite.org/tools/fimo)
"""
pMotifScan                        = Proc(desc = 'Scan motif along the given sequences.')
pMotifScan.input                  = "mfile:file, sfile:file"
pMotifScan.output                 = [
	"outfile:file:{{in.sfile | fn}}-{{args.mids | mid2fn}}.fimo/{{in.sfile | fn}}-{{args.mids | mid2fn}}.bed", 
	"outdir:dir:{{in.sfile | fn}}-{{args.mids | mid2fn}}.fimo"
]
pMotifScan.args.tool              = 'meme'
pMotifScan.args.meme              = params.fimo.value
pMotifScan.args.params            = Box()
pMotifScan.args.mids              = 'all'
pMotifScan.args.pval              = 1e-4
pMotifScan.args.cleanmname        = True
pMotifScan.args.ucsclink          = 'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={}'
pMotifScan.args.nthread           = 1
pMotifScan.tplenvs.mid2fn         = lambda mids: re.sub(r'[^\w_]', '', (''.join(mids) if isinstance(mids, list) else mids))[:16]
pMotifScan.tplenvs.runcmd         = runcmd.py
pMotifScan.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pMotifScan.tplenvs.parallel       = parallel.py
pMotifScan.lang                   = params.python.value
pMotifScan.script                 = "file:scripts/tfbs/pMotifScan.py"


