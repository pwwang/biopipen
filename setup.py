from setuptools import setup, find_packages
import bioprocs

setup (
	name             = 'bioprocs',
	version          = bioprocs.VERSION,
	description      = "A set of procs for bioinformatics using PyPPL",
	url              = "https://github.com/pwwang/bioprocs",
	author           = "pwwang",
	author_email     = "pwwang@pwwang.com",
	license          = "MIT",
	packages         = find_packages()
)
