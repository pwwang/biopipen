import sys
from modkit import Modkit

import cmdy

DEFAULT_CONFIG = dict(
	default = dict(_raise = True),

	bedtools  = dict(_prefix = '-'),
	biobambam = dict(_sep = '=', _prefix = ''),
	bowtie2   = dict(_dupkey = True),
	dtoxog    = dict(_out = cmdy.DEVERR, _prefix = '-'),
	sort      = dict(_sep = '', _dupkey = True),
	gatk      = dict(_dupkey = True),
	liftover  = dict(_prefix = '-', _sep = '='),
	oncotator = dict(_sep = 'auto'),
	maf2vcf   = dict(_sep = ' '),

	# As of picard 2.20.5-SNAPSHOT
	# it's changing in the futher. See: https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
	# Future one should be:
	# picard = dict(_sep = ' ', _prefix = '-')
	picard    = dict(_sep = '=', _prefix = ''),
	plink     = dict(_out = cmdy.DEVERR),
	pyclone   = dict(_raw = True),
	samtools  = dict(_prefix = '-'),
	snpeff    = dict(_prefix = '-'),
	vcfanno   = dict(_prefix = '-'),
)
cmdy.config._load(DEFAULT_CONFIG)

def _modkit_delegate(name):
	return getattr(cmdy, name)

# run command at foreground
fg   = cmdy(_fg = True, _debug = True)
bg   = cmdy(_bg = True, _debug = True)
out  = cmdy(_out = '>')
pipe = cmdy(_pipe = True)

## aliases
rm_rf  = cmdy.rm.bake(r = True, f = True)
ln_s   = cmdy.ln.bake(s = True)
kill_9 = cmdy.kill.bake(s = 9)
wc_l   = cmdy.wc.bake(l = True)

runcmd = lambda cmd: cmdy.bash(c = cmd)

def load_config(conf = None, **kwargs):
	conf = conf or {}
	conf.update(kwargs)
	conf2load = {'default': DEFAULT_CONFIG['default']}
	for key, val in conf.items():
		conf2load[key] = DEFAULT_CONFIG.get(key, {}).copy()
		conf2load[key].update(val if isinstance(val, dict) else {'_exe': val})

	cmdy.config._load(conf2load)
	fg.config._load(conf2load)
	out.config._load(conf2load)

Modkit()
