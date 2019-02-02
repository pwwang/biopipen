from pyppl import Box
from bioprocs.utils import shell

params = {{args.params | repr}}
params.i = {{i.infile | quote}}
params._stdout = {{o.outfile | quote}}

shell.TOOLS.bedtools = {{args.bedtools | quote}}
shell.Shell(subcmd = True, dash = '-', equal = '=').bedtools.merge(**params).run()

