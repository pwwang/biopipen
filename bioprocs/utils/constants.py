import toml
from diot import Diot
from pathlib import Path
from modkit import Modkit

_constantfile = Path(__file__).parent.joinpath('constants.toml')

with _constantfile.open() as fconst:
	_data = Diot(toml.load(fconst))

modkit = Modkit()
modkit.ban('toml', 'Diot', 'Path', 'Modkit', 'modkit', '_constantfile', '_data')
modkit.delegate(lambda name: _data[name])
