import sys
from pyppl import Box
from bioprocs.utils import shell2 as shell

bedtools = {{args.bedtools | quote}}
params   = {{ args.params | repr}}

params.a           = {{i.afile | quote}}
params.b           = {{i.bfiles | repr}}
params.wao         = params.get('wao', True)
params.nonamecheck = params.get('nonamecheck', True)
params._out        = {{o.outfile | quote}}
params._stderr     = sys.stderr

shell.load_config(bedtools = bedtools)

shell.bedtools.intersect(**params)
