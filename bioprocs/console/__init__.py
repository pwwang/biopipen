"""Script entry point for bioprocs"""
import sys
import copy
from pprint import pformat
from diot import Diot
from .. import params
from .utils import highlight_multi, Module, Pipeline, subtract_dict
from .arguments import commands
from pyparam.help import HelpAssembler, Helps
from pyppl.config import config

HELP_ASSEMBLER = HelpAssembler()

def show_params(opts):
    """Show and query a parameter"""
    # get is preserved for Diot
    if opts['get']:
        for pname in opts._:
            if pname not in params:
                raise ValueError('Parameter %r is not defined.' % pname)
            print(params[pname].value)
    else:
        pminfo  = Helps()
        pminfo.add('params', sectype = 'option', prefix = '')
        queries = [query.lower() for query in opts._]
        for key in sorted(params._dict()):
            if key in params._hopts:
                continue
            for query in queries:
                param = params[key]
                if query not in param.str().lower() and query not in key.lower() \
                    and not any(query in desc for desc in param.desc):
                    continue
                pminfo['params'].addParam(param)

        print('\n'.join(highlight_multi(line, queries)
            for line in  HELP_ASSEMBLER.assemble(pminfo)), end = '')

def list_procs(opts):
    """List processes"""
    if not opts._ and not opts.all: # list modules only
        pminfo  = Helps()
        pminfo.add('modules', sectype = 'option', prefix = '')
        for module in Module.modules():
            Module(module).to_helps(pminfo['modules'])
        print('\n'.join(HELP_ASSEMBLER.assemble(pminfo)), end = '')
    elif not opts._:  # list all modules and processes
        pminfo  = Helps()
        modules = Module.modules()
        nmods   = len(modules)
        for module in modules:
            Module(module).to_helps_with_procs(pminfo, nmods)
        print('\n'.join(HELP_ASSEMBLER.assemble(pminfo)), end = '')
    else: # list specific modules and their processes
        pminfo  = Helps()
        modules = Module.modules()
        nmods   = len(opts._)
        for module in opts._:
            if module not in modules:
                raise ValueError('No such module: %s' % module)
            Module(module).to_helps_with_procs(pminfo, nmods)
        print('\n'.join(HELP_ASSEMBLER.assemble(pminfo)), end ='')

def proc(opts):
    """Query processes"""
    modules = Module.modules()
    nmods   = len(modules)
    helps   = Helps()
    for i, module in enumerate(modules):
        module = Module(module)
        if any(q.lower() in module.name.lower() for q in opts._):
            module.to_helps_with_procs(helps, nmods)
        else:
            module.to_helps_as_sec(helps)
            with module.load_procs("%s/%s" % (i+1, nmods)) as procs:
                for name, prc in procs.items():
                    if any(q.lower() in name.lower() for q in opts._):
                        prc.to_helps(helps.select(module.name + ':'))
    print('\n'.join(highlight_multi(line, opts._)
          for line in HELP_ASSEMBLER.assemble(helps)),
          end = '')

def complete(opts):
    """Shell completions"""
    from completions import Completions
    comp = Completions(desc = 'Bioprocs utilities.')
    # add builtin commands
    for cmd, pms in commands._cmds.items():
        if cmd == '_':
            continue
        comp.add_command(cmd, pms._desc)
        pms._add_to_completions(comp.command(cmd), withtype = False, alias = True)

    # add pipelines
    for pipeline in Pipeline.pipelines():
        Pipeline(pipeline).add_to_completions(comp)

    # add modules/processes
    for module in Module.modules():
        Module(module).add_to_completions(comp)

    compcode = comp.generate(shell = opts.s, auto = opts.a)
    if not opts.a:
        print(compcode)

def profile(opts):
    """List avaiable running profiles"""
    if opts._:
        config._load(opts._)

    base = { '_flowchart': {},
        '_log': {},
        'args': {},
        'sgeRunner': {},
        'sshRunner': {},
        'tplenvs': {}}

    helps = Helps()
    helps.add('Profile: "default"', sectype = 'plain')
    helps.select('Profile: "default"').add(
        pformat(subtract_dict(config, base), indent = 2, width = 100))

    default = copy.deepcopy(config)
    for prof in config._profiles:
        if prof == 'default':
            continue
        helps.add('Profile: "%s"' % prof, sectype = 'plain')

        with config._with(prof, copy = True) as profconf:
            profconf = subtract_dict(profconf, default)
            profconf = subtract_dict(profconf, base)
            helps.select('Profile: "%s"' % prof).add(
                pformat(profconf, indent = 2, width = 100))
        config.update(default)

    print('\n'.join(HELP_ASSEMBLER.assemble(helps)), end = '')

def main():
    """Main entry point"""
    pipelines = Pipeline.pipelines()
    command = sys.argv[1] if len(sys.argv) > 1 else None
    if command == 'params':
        _, opts, _ = commands._parse(dict_wrapper = Diot)
        show_params(opts)
    elif command == 'list':
        _, opts, _ = commands._parse(dict_wrapper = Diot)
        list_procs(opts)
    elif command == 'proc':
        _, opts, _ = commands._parse(dict_wrapper = Diot)
        proc(opts)
    elif command in ('completion', 'completions'):
        _, opts, _ = commands._parse(dict_wrapper = Diot)
        complete(opts)
    elif command == 'profile':
        _, opts, _ = commands._parse(dict_wrapper = Diot)
        profile(opts)
    elif command in pipelines:
        # let the pipeline parse the arguments
        Pipeline(command).run()
    elif command in commands._hcmd:
        _, opts, _ = commands._parse(arbi = True, dict_wrapper = Diot)
        if not opts._:
            commands._help(print_and_exit = True)
        if opts._ in commands._cmds:
            commands._cmds[opts._]._help(print_and_exit = True)
        if opts._ in pipelines:
            Pipeline(opts._).module.params._help(print_and_exit = True)
        if '.' in opts._: # process
            module, prc = opts._.split('.')
            try:
                Module(module).procs()[prc].printHelps()
            except KeyError:
                raise KeyError('Module %r does not have process: %s' % (module, prc)) from None
        else: # assume module
            list_procs(Diot(_ = [opts._]))
    elif not sys.argv[0]:
        raise RuntimeError('This package has to run as a command line tool.')
    elif not command:
        commands._help(print_and_exit = True)
    else:
        # Hold helps first, because we haven't assembled the help page yet
        hopts = commands[command]._hopts
        commands[command]._hopts = []
        command, opts, _ = commands._parse(arbi = True, dict_wrapper = Diot)
        # resume it
        commands[command]._hopts = hopts
        if '.' in command:
            module, prc = command.split('.')
            module = Module(module)
            if prc not in module.procs():
                raise ValueError('No such process: %s' % command)
            module.procs()[prc].run(opts)
        else: # listing module processes
            list_procs(Diot(_ = [command]))
