from setuptools import setup, find_packages
import bioprocs

setup (
	name             = 'bioprocs',
	version          = bioprocs.VERSION,
	description      = "A set of procs for bioinformatics using pyppl",
	url              = "https://github.com/pwwang/pyppl",
	author           = "pwwang",
	author_email     = "pwwang@pwwang.com",
	license          = "MIT",
	packages         = find_packages(),
	classifiers      = [
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: BSD",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.2",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
    ]
)