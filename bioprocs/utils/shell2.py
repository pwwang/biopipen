import sys
from modkit import Modkit

import cmdy

cmdy.config._load(dict(
	default = dict(_raise = False),
	sort    = dict(_sep = '', _dupkey = True),

	# As of picard 2.18.27-SNAPSHOT
	# it's changing in the futher. See: https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
	# Future one should be:
	# picard = dict(_sep = ' ', _prefix = '-')
	picard    = dict(_sep = '=', _prefix = ''),
	plink     = dict(_out = cmdy.DEVERR),
	oncotator = dict(_sep = 'auto'),
	bowtie2   = dict(_dupkey = True),
	gatk      = dict(_dupkey = True),
	biobambam = dict(_sep = '=', _prefix = ''),
	bedtools  = dict(_prefix = '-'),
	samtools  = dict(_prefix = '-'),
))


def _modkit_delegate(name):
	return getattr(cmdy, name)

# run command at foreground
fg   = cmdy(_fg = True, _raise = True)
out  = cmdy(_out = '>')

## aliases
rm_rf  = cmdy.rm.bake(r = True, f = True)
ln_s   = cmdy.ln.bake(s = True)
kill_9 = cmdy.kill.bake(s = 9)

runcmd = lambda cmd: cmdy.bash(c = cmd)

def load_config(conf = None, **kwargs):
	conf = conf or {}
	conf.update(kwargs)
	conf2load = {
		key: (val if isinstance(val, dict) else {'_exe': val})
		for key, val in conf.items()
	}
	cmdy.config._load(conf2load)
	fg.config._load(conf2load)
	out.config._load(conf2load)

Modkit()