#!/usr/bin/env python
# Using PyPPL to distribute and run commands.
# @params:
#	`cmds`  : The cmd list.
#	`runner`: The runner.
#	`forks` : How many jobs to run simultaneously.
import sys
from pyppl import PyPPL, Proc
from bioprocs import params

params.prefix('-')
params.cmds.desc   = 'The cmd list. If not provided, STDIN will be used.'
params.runner      = 'local'
params.runner.desc = 'The runner.'
params.forks       = 1
params.forks.desc  = 'How many jobs to run simultaneously.'


def _proc(kwargs = None):
	global params

	kwargs = kwargs or {}
	if kwargs and '' in params._props['hopts']:
		del params._props['hopts'][params._props['hopts'].index('')]
		
	for key, val in kwargs.items():
		setattr(params, key, val)

	params = params.parse(raiseExc = bool(kwargs), args = [] if kwargs else sys.argv).asDict()
	if not params['cmds']:
		raise ValueError('Option -cmds is required.')

	if not isinstance(params['cmds'], list):
		params['cmds'] = params['cmds'].splitlines()

	pCmdRunner = Proc(desc = 'Using PyPPL to distribute and run commands.')
	pCmdRunner.runner = params.get('runner', 'local')
	pCmdRunner.cache  = False
	pCmdRunner.forks  = int(params.get('forks', 1))
	pCmdRunner.input  = {'cmd': params['cmds']}
	pCmdRunner.script = '{{in.cmd}}'
	return pCmdRunner

def run(*args, **kwargs):
	proc = _proc(kwargs)
	PyPPL({'log':{'file':None}}).start(proc).run()

if __name__ == '__main__':
	run()
