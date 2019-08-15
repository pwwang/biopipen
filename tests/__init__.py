from pathlib import Path
from remotedata import remotedata
remotedata.manager.cachedir   = Path(__file__).parent / 'testdata'
remotedata.manager.conf.repos = 'pwwang/bioprocs-testdata'

def assertInfile(file, *strings):
	content = Path(file).read_text()
	for string in strings:
		assert string in content
		content = content[content.find(string) + len(string):]
