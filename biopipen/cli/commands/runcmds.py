#!/usr/bin/env python
"""Using PyPPL to distribute and run commands."""
import sys
from os import path
from tempfile import gettempdir
from pyppl import PyPPL, Proc
from diot import Diot
from pyparam import Params
from biopipen import params as bp_params

help_group = 'SCRIPTS'

params = Params(desc=__doc__, help_on_void=False)

params.add_param(
    'cmds',
    type=str,
    type_frozen=False,
    desc=('The cmd list. If not provided, STDIN will be used.'
          'To enter list of cmd, one should overwrite the type by --cmds:list')
)
params.add_param('runner', default='local', desc='The runner', force=True)
params.add_param('intype', type='choice', choices=['stdin', 'cmds', 'file'],
                 force=True,
                 desc='Type of `cmds`. One of {choices}')
params.add_param('cache', default=False, force=True,
                 desc='Cache the jobs or not?')
params.add_param('ppldir', default='./workdir', type='path', force=True,
                 desc='The ppldir to save the pipeline data.')
params.add_param('forks', default=1, force=True,
                 desc='How many jobs to run simutaneously.')

def main(opts):
    """Main function"""
    opts = bp_params.parse() | opts

    if not isinstance(opts.cmds, list):
        if opts.intype == 'cmds':
            opts.cmds = opts.cmds.splitlines()
        elif opts.intype == 'file':
            with open(opts.cmds) as fcmd:
                opts.cmds = fcmd.read().splitlines()
        else:
            opts.cmds = []
            for line in sys.stdin:
                opts.cmds.append(line.strip())
    # pylint: disable=invalid-name
    pCmdRunner = Proc(desc='Using PyPPL to distribute and run commands.')
    pCmdRunner.runner = opts.get('runner', 'local')
    pCmdRunner.cache = opts.cache
    pCmdRunner.errhow = 'halt'
    pCmdRunner.forks = int(opts.get('forks', 1))
    pCmdRunner.input = {'cmd': opts.cmds}
    pCmdRunner.output = 'output:var:0'
    pCmdRunner.script = '{{i.cmd}}'

    PyPPL(ppldir=opts.ppldir).start(pCmdRunner).run()


if __name__ == "__main__":
    main(params.parse())
