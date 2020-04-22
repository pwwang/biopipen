"""A set of bioinformatics for PyPPL"""

__version__ = '0.1.2'

import inspect
from pathlib import Path
from tempfile import gettempdir
from sys import executable, modules
from pyppl import Proc
from pyparam import Params
from varname import varname
from . import _envs

# pylint: disable=invalid-name
# open to R (reticulate) to get the path of r util scripts
HERE = Path(__file__).resolve().parent

params = Params()
params._load_file(HERE / 'params.toml')

params.cachedir.value = Path(params.cachedir.value).expanduser().as_posix()
params.tmpdir.value = (gettempdir()
                       if Path(params.tmpdir.value) == Path('/tmp')
                       else params.tmpdir.value)

cfgfiles = [
    Path.home() / '.bioprocs.config',  # values overwritten
    Path.home() / '.bioprocs.json',
    Path.home() / '.bioprocs.toml',
    Path('.') / '.bioprocs.config',
    Path('.') / '.bioprocs.json',
    Path('.') / '.bioprocs.toml',
    'bioprocs.osenv'
]
for cfgfile in cfgfiles:
    if isinstance(cfgfile, Path) and not cfgfile.exists():
        continue
    params._load_file(cfgfile)

# lock the params in case the options are overwritten unintentionally.
if not params._locked:
    params._locked = True

cachedir = Path(params.cachedir.value)
if not cachedir.exists():
    cachedir.mkdir()

EXT_MAP = {
    'Rscript': 'R',
    'python': 'py',
    'python2': 'py',
    'python3': 'py',
}


def _findscript(script, callerdir):
    if not script or not script.startswith('file:'):
        return script
    scriptfile = Path(script[5:])
    if scriptfile.is_absolute():
        return script
    scriptfile = callerdir.joinpath(scriptfile)
    return 'file:{}'.format(scriptfile)

# pylint: disable=redefined-builtin,invalid-name
def proc_factory(id=None, tag='notag', desc='No description.', **kwargs):
    """A factory to produce processes with default script,
    envs and report_template"""
    # in case if we have too long description
    id = id or varname()
    proc = Proc(id, tag=tag, desc=desc, **kwargs)
    lang = Path(proc.lang).name
    ext = '.' + EXT_MAP.get(lang, lang)
    callerfile = Path(
        inspect.getframeinfo(inspect.currentframe().f_back).filename)

    script = proc._script or 'file:scripts/{}/{}{}'.format(
        callerfile.stem, id, ext)
    proc.script = _findscript(script, callerfile.parent)
    report_template = proc.config.report_template or \
     'file:reports/{}/{}.md'.format(callerfile.stem, id)
    report_template = _findscript(report_template, callerfile.parent)
    if Path(report_template[5:]).is_file() or (report_template and
                                               report_template[:5] != 'file:'):
        proc.config.report_template = report_template
    return proc
