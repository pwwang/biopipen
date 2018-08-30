from setuptools import setup, find_packages
# get version
from os import path
verfile = path.join(path.dirname(__file__), 'bioprocs', '__init__.py')
with open(verfile) as vf:
	VERSION = vf.readline().split('=')[1].strip()[1:-1]

setup (
	name             = 'bioprocs',
	version          = VERSION,
	description      = "A set of procs for bioinformatics using PyPPL",
	url              = "https://github.com/pwwang/bioprocs",
	author           = "pwwang",
	author_email     = "pwwang@pwwang.com",
	license          = "Apache License Version 2.0",
	packages         = find_packages(),
	scripts          = ['bin/bioprocs'],
	install_requires = [
		'PyPPL'
	],
	include_package_data = True,
	classifiers          = [
		"Intended Audience :: Developers",
		"Natural Language :: English",
		"Operating System :: MacOS :: MacOS X",
		"Operating System :: POSIX",
		"Operating System :: POSIX :: BSD",
		"Operating System :: POSIX :: Linux",
		"Programming Language :: Python"
	]
	
)
