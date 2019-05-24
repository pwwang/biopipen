from sys import stderr
from pyppl import Box
from bioprocs.utils import shell2 as shell

params = {{args.params | repr}}

params['a']           = {{i.afile | quote}}
params['b']           = {{i.bfile | quote}}
params['wa']          = params.get('wa', True)
params['wb']          = params.get('wb', True)
params['nonamecheck'] = params.get('nonamecheck', True)
params['_out']        = {{o.outfile | quote}}
params['_stderr']     = stderr

shell.load_config(bedtools = {{args.bedtools | quote}})
shell.bedtools.intersect(**params)
