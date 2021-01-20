"""Initialize params"""
from os import environ
from pathlib import Path
from tempfile import gettempdir

import modkit
from pipen_args import args as params
from pyparam import Namespace

from .defaults import PKG_NAME

# pylint: disable=invalid-name
# open to R (reticulate) to get the path of r util scripts
HERE = Path(__file__).parent.resolve()

params.from_file(HERE / 'args.toml', show=False)

cfgfiles = [
    Path.home() / f'.{PKG_NAME}.toml',
    Path('.') / f'.{PKG_NAME}.toml',
    f'{PKG_NAME}.osenv'
]
for cfgfile in cfgfiles:
    if isinstance(cfgfile, Path) and not cfgfile.exists():
        continue
    params.from_file(cfgfile, force=True, show=False)

# expanduser if cachedir has ~
params.get_param('cachedir').callback = lambda val: str(Path(val).expanduser())

# make sure we are using system's tmpdir, rather than /tmp
# sometimes when peope say /tmp, they are actually mean system's tmpdir
# set by $TMPDIR for example.
params.get_param('tmpdir').callback = lambda val: (
    val if val
    else environ['TMPDIR'] if 'TMPDIR' in environ
    else gettempdir()
)

args = Namespace()
params.values(args)

Path(args.cachedir).mkdir(exist_ok=True)
Path(args.tmpdir).mkdir(exist_ok=True)

def parse(cli_args=None):
    params.parse(cli_args)
    params.values(args)
    return args

def __getattr__(name):
    if name == 'parsed':
        return parse([])
    raise AttributeError(name)

modkit.install(__name__)
