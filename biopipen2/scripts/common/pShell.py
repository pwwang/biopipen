from diot import Diot
from biopipen.utils import ensureDiot, shell2 as shell

cmd = {{args.cmd | repr}}
if not cmd:
	raise ValueError('No cmd (args.cmd) specified')

cmd = {{args.cmd | render | quote}}

shell.bash(c = cmd).fg
