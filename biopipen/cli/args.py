"""Arguments of bioprocs"""
import logging
from pathlib import Path
from pyparam import Params
from .modules import Module
from .commands import Command
from .utils import file_newer, logger, ROOT_MODULE

CACHED_ARGS = Path(f'~/.{ROOT_MODULE}.args.toml').expanduser()

def help_modifier(help_param, help_cmd): # pylint: disable=unused-argument
    """Modify the help argument or command"""
    help_cmd.group = 'Other commands'

def help_callback(helps):
    """Modify help page"""
    # don't show all processes
    del helps['MODULES'][:]
    helps['MODULES'].append((['{modname}'], ['List process of the module.']))

params = Params(prog=ROOT_MODULE,
                desc='A set of procs for bioinformatics.',
                help_callback=help_callback,
                help_modifier=help_modifier)
# see if anything changed, otherwise load from the CACHED_ARGS
if (not CACHED_ARGS.is_file() or
        not Command.use_cache(CACHED_ARGS) or
        not Module.use_cache(CACHED_ARGS) or
        file_newer(__file__, CACHED_ARGS)):
    logger.setLevel(logging.DEBUG)
    logger.debug('# Loading and caching parameters ...')
    Command.to_params(params)
    Module.to_params(params)
    logger.setLevel(logging.INFO)
    params.to_file(CACHED_ARGS)
else:
    params.from_file(CACHED_ARGS)
