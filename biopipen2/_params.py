"""Initialize params"""
from pathlib import Path
from tempfile import gettempdir
from pyparam import Params, Namespace as PMNamespace

class Namespace(PMNamespace):
    """Give default value '' when name not exists"""

    def __getattr__(self, name):
        try:
            return super().__getattr__(name)
        except AttributeError:
            from .cli.utils import logger
            logger.warning("Parameter %r not defined, using ''.", name)
            return ''

# pylint: disable=invalid-name
# open to R (reticulate) to get the path of r util scripts
HERE = Path(__file__).resolve().parent
ROOT_MODULE = 'biopipen'

params = Params(help_on_void=False)
params.from_file(HERE / '_params.toml', show=False)

cfgfiles = [
    Path.home() / f'.{ROOT_MODULE}.config',  # values overwritten
    Path.home() / f'.{ROOT_MODULE}.json',
    Path.home() / f'.{ROOT_MODULE}.toml',
    Path('.') / f'.{ROOT_MODULE}.config',
    Path('.') / f'.{ROOT_MODULE}.json',
    Path('.') / f'.{ROOT_MODULE}.toml',
    f'{ROOT_MODULE}.osenv'
]
for cfgfile in cfgfiles:
    if isinstance(cfgfile, Path) and not cfgfile.exists():
        continue
    params.from_file(cfgfile, force=True, show=False)

# expanduser if cachedir has ~
param_cachedir = params.get_param('cachedir')
if param_cachedir:
    param_cachedir.callback = lambda val: str(Path(val).expanduser())

# make sure we are using system's tmpdir, rather than /tmp
# sometimes when peope say /tmp, they are actually mean system's tmpdir
# set by $TMPDIR for example.
param_tmpdir = params.get_param('tmpdir')
if param_tmpdir:
    param_tmpdir.callback = (
        lambda val: gettempdir() if val in ('/tmp', '/tmp/') else val
    )

opts = Namespace()
params.values(opts)
if not Path(opts.cachedir).exists():
    Path(opts.cachedir).mkdir()

def update_params():
    """Update opts when paramters modified in params"""
    params.value(opts)
