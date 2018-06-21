#!/usr/bin/env python
# Using PyPPL to distribute and run commands.
import sys
from os import path
from pyppl import PyPPL, Proc
from bioprocs import params

params.prefix('-')
params.cmds.required = True
params.cmds.desc     = 'The cmd list. If not provided, STDIN will be used.'
params.runner        = 'local'
params.runner.desc   = 'The runner.'
params.forks         = 1
params.forks.desc    = 'How many jobs to run simultaneously.'

def _procconfig(kwargs = None):
	params = _getparams(kwargs or {})

	if not isinstance(params['cmds'], list):
		params['cmds'] = params['cmds'].splitlines()

	pCmdRunner = Proc(desc = 'Using PyPPL to distribute and run commands.')
	pCmdRunner.runner = params.get('runner', 'local')
	pCmdRunner.cache  = False
	pCmdRunner.forks  = int(params.get('forks', 1))
	pCmdRunner.input  = {'cmd': params['cmds']}
	pCmdRunner.script = '{{in.cmd}}'

	config = {'proc': {'file': None}}
	return pCmdRunner, config

def _getparams(kwargs):
	global params
	
	if len(sys.argv) > 1 and sys.argv[1] == path.splitext(path.basename(__file__))[0]:
		# called from api
		if '' in params._props['hopts']:
			del params._props['hopts'][params._props['hopts'].index('')]
		for key, val in kwargs.items():
			setattr(params, key, val)
		return params.parse(args = []).asDict()
	else:
		# called directly
		return params.parse().asDict()

def run(*args, **kwargs):
	proc, config = _procconfig(kwargs)
	PyPPL(config).start(proc).run()

if __name__ == '__main__':
	run()
