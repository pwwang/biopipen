"""Generate completion code to be integrated into the shell.
Completion code should be saved (appended) to:
* bash: `~/.bashrc`
* zsh `~/.zshrc`
* fish `~/.config/fish/completions/<script>.fish`
"""
import sys
from pathlib import Path
from pyparam import Params
from ..modules import MODULE_FILES, MODULE_DIR
from ..proxy import Module
from ..utils import logger, ROOT_MODULE

params = Params(desc=__doc__)
params.add_param('shell', type='choice', choices=['bash', 'fish', 'zsh'],
                 desc='The shell name. One of {choices}', required=True)
params.add_param('module', type=str, default='self',
                 desc=['Which module to generate completion code for.',
                       f'* `self`: generate for `{ROOT_MODULE}`',
                       f'* `all`: generate for all modules, including `self`',
                       '  [yellow]WARNING[/yellow]: codes will be '
                       'automatically appended ',
                       '  to target file(s).',
                       f'* `{{{{module}}}}`: generate for `{{{{module}}}}`',
                       '  For fish, code should be saved to ',
                       '  `~/.config/fish/completions/'
                       f'{ROOT_MODULE}-{{{{module}}}}.fish`'])

def _saving(code, shell, name=None):
    """Saving code to corresponding file"""
    def saving_to(code, path):
        print(code)
        logger.warning('Saving code to %r', str(path))
        with open(path.expanduser(), 'a') as file:
            file.write('\n' + code + '\n')

    if shell == 'zsh':
        saving_to(code, Path('~/.zshrc'))
    elif shell == 'bash':
        saving_to(code, Path('~/.bashrc'))
    else:
        saving_to(code, Path(f'~/.config/fish/completions/{name}.fish'))

def main(opts):
    """Main function"""
    if opts.module == 'self':
        from ..args import params as root_params
        print(root_params.shellcode(opts.shell))
    elif opts.module == 'all':
        from ..args import params as root_params
        code = root_params.shellcode(opts.shell)
        _saving(code, opts.shell, ROOT_MODULE)
        for file in sorted(MODULE_FILES):
            module = Module(file)
            _saving(module.params.shellcode(opts.shell), opts.shell,
                    f'{ROOT_MODULE}-{module.name}')

    elif opts.module not in (file.stem for file in MODULE_FILES):
        logger.error("No such module, use `%s list` to show available modules",
                     ROOT_MODULE)
        sys.exit(1)
    else:
        module = Module(MODULE_DIR.joinpath(f"{opts.module}.py"))
        print(module.params.shellcode(opts.shell))

if __name__ == '__main__':
    main(params.parse())
