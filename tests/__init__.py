from pathlib import Path
from remotedata import remotedata
remotedata.manager.cachedir   = Path(__file__).parent / 'testdata'
remotedata.manager.conf.repos = 'pwwang/bioprocs-testdata'
