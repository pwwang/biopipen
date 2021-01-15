import toml
from diot import Diot
from pathlib import Path
from modkit import modkit
modkit.ban('toml', 'Diot', 'Path', 'Modkit', 'modkit', '_constantfile', '_data')

_constantfile = Path(__file__).parent.joinpath('constants.toml')

with _constantfile.open() as fconst:
	_data = Diot(toml.load(fconst))

@modkit.delegate
def _constant_delegate(module, name):
	if name in _data:
		return _data[name]
	raise AttributeError
