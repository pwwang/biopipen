from pyppl import Diot
from bioprocs.utils import ensureDiot, shell2 as shell

cmd = {{args.cmd | repr}}
if not cmd:
	raise ValueError('No cmd (args.cmd) specified')

cmd = {{args.cmd | render | quote}}

shell.fg.bash(c = cmd)
