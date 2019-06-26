from pyppl import Box
from bioprocs.utils import ensureBox, shell2 as shell

args = {{i.args if isinstance(i.args, str) else repr(i.args)}}
args = ensureBox(args)

cmd = {{args.cmd | repr}}

if not cmd:
	raise ValueError('No cmd (args.cmd) specified')

if isinstance(cmd, dict):
	cmd, exe = cmd.items()[0]
elif isinstance(cmd, (tuple, list)):
	cmd, exe = cmd
else:
	cmd, exe = cmd, cmd

cmds    = cmd.split('.')
command = shell.fg
for cmd in cmds:
	command = getattr(command, cmd)
command(_exe = exe, **args)
