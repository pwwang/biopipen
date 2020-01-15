import pytest
from pathlib import Path

def assertInfile(file, *strings):
	content = Path(file).read_text()
	for string in strings:
		assert string in content
		content = content[content.find(string) + len(string):]
