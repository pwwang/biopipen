"""Add a module to pyproject.toml as a console script. Internal use ONLY."""
import sys
from pathlib import Path
from pyparam import Params, POSITIONAL
from ..modules import MODULE_FILES
from ..utils import logger, ROOT_MODULE

params = Params(desc=__doc__)
params.add_param(POSITIONAL, type=str, required=True,
                 desc='Name of the module to add.')

PYPROJECT_TOML = Path(__file__).parent.parent.parent.parent.joinpath(
    'pyproject.toml'
)

def _scan_scripts():
    # None: nothing hit
    # start: hit [tool.poetry.scripts]
    # end: when status is start and hit empty line or [
    status = None
    scripts = []
    before = []
    after = []
    with open(PYPROJECT_TOML) as file:
        for line in file:
            line = line.strip()
            if not status:
                before.append(line)
                if line == '[tool.poetry.scripts]':
                    status = 'start'

            elif status == 'start' and line.startswith(ROOT_MODULE + '-'):
                scripts.append(line.split('=', 1)[0].split('-')[1].strip())

            elif status == 'start' and (not line or line[:1] == '['):
                after.append(line)
                status = 'end'

            elif status == 'start':
                before.append(line)

            elif status == 'end':
                after.append(line)

    return before, scripts, after

_SCANNED_SCRIPTS = _scan_scripts()
BEFORE_SCRIPTS = _SCANNED_SCRIPTS[0]
SCRIPTS = _SCANNED_SCRIPTS[1]
AFTER_SCRIPTS = _SCANNED_SCRIPTS[2]

def _module_to_add(modname):
    if modname in SCRIPTS:
        logger.warning('Module has already been added: %r', modname)
    else:
        if modname not in (file.stem for file in MODULE_FILES):
            logger.error('Module does not exist: %r', modname)
        else:
            ret = (f"{ROOT_MODULE}-{modname} = "
                   f"'{ROOT_MODULE}.cli.proxy:{modname}'")
            logger.info('Adding %r', modname)
            return ret
    return None

def main(opts):
    """Main entry point"""
    if not PYPROJECT_TOML.is_file():
        logger.error('This command is for internal use ONLY.')
        sys.exit(1)

    modname = opts[POSITIONAL]
    lines = []
    if modname == 'all':
        for file in MODULE_FILES:
            modname = file.stem
            line = _module_to_add(modname)
            if line:
                lines.append(line)
    if lines:
        Path(str(PYPROJECT_TOML) + '.bak').write_text(
            PYPROJECT_TOML.read_text()
        )
        with open(PYPROJECT_TOML, 'w') as file:
            file.write('\n'.join(BEFORE_SCRIPTS) + '\n')
            file.write('\n'.join(lines) + '\n')
            file.write('\n'.join(AFTER_SCRIPTS) + '\n')

if __name__ == '__main__':
    main(params.parse())
