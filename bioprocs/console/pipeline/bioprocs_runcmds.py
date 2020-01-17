#!/usr/bin/env python
"""Using PyPPL to distribute and run commands."""
import sys
from os import path
from tempfile import gettempdir
from pyppl import PyPPL, Proc
from diot import Diot
from bioprocs import params

params._prefix = '-'
params._hbald = False
params.cmds.desc = 'The cmd list. If not provided, STDIN will be used.'
params.runner = 'local'
params.runner.desc = 'The runner.'
params.intype = 'stdin'  # or file, stdin, or cmds (pass cmds directly)
params.intype.desc = 'Type of option `cmds`, cmds or file?'
params.cache = False
params.cache.desc = 'Cache the jobs or not?'
params.ppldir = path.join(gettempdir(), 'bioprocs.workdir')
params.ppldir.desc = 'The ppldir to save the pipeline data.'
params.forks = 1
params.forks.desc = 'How many jobs to run simultaneously.'


def main():
    """Main function"""
    opts = params._parse(dict_wrapper=Diot)

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
    pCmdRunner.script = '{{i.cmd}}'

    PyPPL(ppldir=opts.ppldir).start(pCmdRunner).run()


if __name__ == "__main__":
    main()
