import warnings

from pipen import Proc as PipenProc
from .envs import envs
from .args import params
from .utils import logger

# Warn if params are not parsed
if not any(key in params.params for key in params.help_keys):
    logger.warning('Params not parsed, using all default values.')
    logger.warning('To suppress this warnings, try '
                   '`from biopipen.core.args import parsed`')

class Proc(PipenProc):

    envs = envs
