import pytest
from pathlib import Path
from remotedata import remotedata
from simpleconf import Config

@pytest.fixture
def config():
	ret = Config()
	ret._load({
		'default': {
			'cachedir': Path(__file__).resolve().parent / 'testdata',
			'source': 'github',
			'repos': 'pwwang/bioprocs-testdata'
		}
	}, '~/.remotedata.yaml', 'REMOTEDATA.osenv')
	return ret

@pytest.fixture
def rdata(config):
	return remotedata(config)