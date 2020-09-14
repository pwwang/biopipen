"""Load modules"""
from typing import Union
from functools import lru_cache
from pathlib import Path
from rich import print as rich_print, box
from rich.table import Table
from .utils import (
    file_newer, CommandABC, logger
)

MODULE_DIR = Path(__file__).parent.parent
MODULE_FILES = [path for path in MODULE_DIR.glob('*.py')
                if not path.name.startswith('_')]

def _read_doc(path):
    """Read the doc of a module file"""
    start_tag = None
    doc = []
    with open(path) as file:
        for line in file:
            line = line.strip()
            if start_tag and not line:
                # hit empty line, use the first paragraph
                break
            if start_tag and line.endswith(start_tag):
                if line[:-3]:
                    doc.append(line[:-3])
                break
            if start_tag:
                doc.append(line)
            elif line[:3] not in ('"""', "'''"):
                continue
            else: # not start_tag, line[:3] in (...)
                if line[:3] == line[-3:] and len(line) > 6:
                    doc.append(line[3:-3])
                    break
                start_tag = line[:3]
                if line[3:]:
                    doc.append(line[3:])
    return "\n".join(doc)

class Module(CommandABC):
    """A module of biopipen"""

    @classmethod
    @lru_cache()
    def collections(cls):
        """Load all modules"""
        ret = {}
        logger.debug('- Loading sub-modules ...')
        for modfile in MODULE_FILES:
            logger.debug('  Loading %s', modfile.stem)
            ret[modfile.stem] = cls(modfile)
        return ret

    @classmethod
    def use_cache(cls, cache_file: Path) -> bool:
        """Tell should we use the cache file"""
        return file_newer(cache_file, MODULE_FILES)

    @classmethod
    def is_a(cls, name: str) -> bool:
        """Tell if a given name is a module"""
        return MODULE_DIR.joinpath(name + '.py').is_file()

    @classmethod
    def get(cls, name: str) -> Union["Module", "Process"]:
        """Get the module or process instance by name"""
        return cls(MODULE_DIR.joinpath(name + '.py'))

    @classmethod
    def to_params(cls, params: "Params"):
        """Load the modules to params"""
        for modname, module in cls.collections().items():
            params.add_command(modname,
                               desc=module.doc,
                               group='MODULES',
                               help_on_void=False)

    def __init__(self, path: Path):
        self.path = path
        self.name = path.stem
        # load the module from path
        self.doc = _read_doc(path) or '[ Undocumented ]'

    def run(self, opts):
        """Nothing to do but just list all processes"""
        from .proxy import Module as ProxyModule
        module = ProxyModule(self.path)

        table = Table(box=box.SIMPLE)
        table.add_column('Process', style="green")
        table.add_column('Description')

        for name, proc in module.procs.items():
            table.add_row(name,
                          (proc.doc.description or proc.proc.desc).rstrip())

        rich_print(table)
