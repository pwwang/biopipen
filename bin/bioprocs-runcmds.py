#!/usr/bin/env python
# Using PyPPL to distribute and run commands.
import sys
from os import path
from pyppl import PyPPL, Proc
from bioprocs import params
from tempfile import gettempdir

#params.prefix('-')
params.cmds.required = True
params.cmds.desc     = 'The cmd list. If not provided, STDIN will be used.'
params.runner        = 'local'
params.runner.desc   = 'The runner.'
params.intype        = 'stdin' # or file, stdin, or cmds (pass cmds directly)
params.intype.desc   = 'Type of option `cmds`, cmds or file?'
params.cache         = False
params.cache.desc    = 'Cache the jobs or not?'
params.ppldir        = path.join(gettempdir(), 'bioprocs.workdir')
params.ppldir.desc   = 'The ppldir to save the pipeline data.'
params.forks         = 1
params.forks.desc    = 'How many jobs to run simultaneously.'

params = params.parse()

if not isinstance(params.cmds, list):
	if params.intype == 'cmds':
		params.cmds = params.cmds.splitlines()
	elif params.intype == 'file':
		with open(params.cmds) as f:
			params.cmds = f.read().splitlines()
	else:
		params.cmds = []
		for line in sys.stdin:
			params.cmds.append(line.strip())

pCmdRunner        = Proc(desc = 'Using PyPPL to distribute and run commands.')
pCmdRunner.runner = params.get('runner', 'local')
pCmdRunner.cache  = params.cache
pCmdRunner.errhow = 'halt'
pCmdRunner.forks  = int(params.get('forks', 1))
pCmdRunner.input  = {'cmd': params.cmds}
pCmdRunner.script = '{{i.cmd}}'

config = {'_log': {'file': None}, 'default': {'ppldir': params.ppldir}}
PyPPL(config).start(pCmdRunner).run()
