VERSION = "0.0.1a"

import json
from os import path
from sys import modules, stderr
from pyppl import params

DEFAULTS = {
	# tools
	"curl"     : "curl",
	"curl.desc": "The path of command line tool curl.",

	# langs
	"python"      : "python",
	"python.desc" : "The path of python.",
	"Rscript"     : "Rscript",
	"Rscript.desc": "The path of Rscript.",
	"bash"        : "bash",
	"bash.desc"   : "The path of bash.",
}

params.loadDict(DEFAULTS)

cfgfiles = [
	path.join (path.expanduser('~'), ".bioProcs"),   # values overwritten
	path.join (path.expanduser('~'), ".bioProcs.json")
]
for cfgfile in cfgfiles:
	if not path.exists(cfgfile):
		continue
	params.loadFile (cfgfile)