"""Manage command for biopipen"""
from typing import Union
from functools import lru_cache
from pathlib import Path
from ..utils import (
    ROOT_MODULE,
    file_newer,
    CommandABC,
    logger,
    import_from_path
)

COMMAND_DIR = Path(__file__).parent
COMMAND_FILES = [path for path in COMMAND_DIR.glob('*.py')
                 if not path.name.startswith('_')]

class Command(CommandABC):
    """Manage commands"""

    @classmethod
    @lru_cache()
    def collections(cls):
        """Load all commands"""
        ret = {}
        logger.debug('- Loading commands ...')
        for cmdfile in COMMAND_FILES:
            logger.debug('  Loading %s', cmdfile.stem)
            ret[cmdfile.stem] = cls(cmdfile)
        return ret

    @classmethod
    def use_cache(cls, cache_file: Path) -> bool:
        """Tell should we use the cache file"""
        return file_newer(cache_file, COMMAND_FILES)

    @classmethod
    def is_a(cls, name: str) -> bool:
        return COMMAND_DIR.joinpath(name + '.py').is_file()

    @classmethod
    def get(cls, name: str) -> "Command":
        return cls(COMMAND_DIR.joinpath(name + '.py'))

    @classmethod
    def to_params(cls, params):
        from rich import print
        for cmdname, command in cls.collections().items():
            command.params.names = [cmdname]
            params.add_command(command.params, group=command.help_group)

    def __init__(self, cmdfile: Path):
        self.module = import_from_path(cmdfile)

    def run(self, opts):
        self.module.main(opts)

    @property
    def params(self):
        """Get params of the command script"""
        return self.module.params

    @property
    def help_group(self):
        """Get the help group if defined in the script"""
        try:
            return self.module.help_group
        except AttributeError:
            return 'COMMANDS'
