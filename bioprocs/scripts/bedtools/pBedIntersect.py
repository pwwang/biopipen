from pyppl import Box
from bioprocs.utils import shell

params = {{args.params | repr}}

params['a']           = {{i.afile | quote}}
params['b']           = {{i.bfile | quote}}
params['wao']         = params.get('wao', True)
params['nonamecheck'] = params.get('nonamecheck', True)
params['_stdout']     = {{o.outfile | quote}}

shell.TOOLS.bedtools = {{ args.bedtools | quote}}
shell.Shell(subcmd = True, dash = '-', equal = ' ').bedtools.intersect(**params).run()

