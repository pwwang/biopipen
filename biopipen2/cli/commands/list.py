"""List all modules or processes of a given module"""
import sys
from pyparam import Params, POSITIONAL
from rich import print as rich_print, box
from rich.table import Table
from ..utils import logger

params = Params(desc=__doc__, help_on_void=False)
params.add_param(POSITIONAL, type=str,
                 desc='List processes of the module. If not given, '
                 'will list all modules.')

def main(opts):
    from ..args import params
    if not opts[POSITIONAL]:
        table = Table(box=box.SIMPLE)
        table.add_column('Module', style="green")
        table.add_column('Description')

        # for modname, module in Module.collections().items():
        #     table.add_row(modname, module.doc.rstrip())

        # use the cache instead
        for module_command in params.command_groups['MODULES']:
            table.add_row(module_command.names[0],
                          module_command.desc[0].rstrip())

        rich_print(table)
    else:
        if opts[POSITIONAL] not in (
                cmd.names[0] for cmd in params.command_groups['MODULES']
        ):
            logger.error('No such module %r. Available modules:',
                         opts[POSITIONAL])
            main({POSITIONAL: None})
            sys.exit(1)

        table = Table(box=box.SIMPLE)
        table.add_column('Process', style="green")
        table.add_column('Description')
        for proc in params.command_groups['PROCESSES']:
            procname = proc.names[0]
            if not procname.startswith(opts[POSITIONAL] + '.'):
                continue
            procname = procname[len(opts[POSITIONAL]) + 1:]

            table.add_row(procname, proc.desc[0].rstrip())

        rich_print(table)

if __name__ == '__main__':
    parsed = params.parse()
    main(parsed)
