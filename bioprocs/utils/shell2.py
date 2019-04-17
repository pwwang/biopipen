# requires https://github.com/pwwang/sh
# install: pip install git+https://github.com/pwwang/sh
from pyppl import Box
from sh import BAKED_ARGS
# defaults
BAKED_ARGS.update(dict(
	sort = dict(_sep = '', _duplistkey = True),

	# As of picard 2.18.27-SNAPSHOT
	# it's changing in the futher. See: https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
	# Future one should be:
	# picard = dict(_sep = ' ', _prefix = '-')
	picard = dict(_sep = '=', _prefix = ''),
	oncotator = dict(_sep = 'auto'),
))

def update_args(*args, **kwargs):
	for arg in args:
		assert isinstance(arg, dict)
		update_args(**arg)
	for k, v in kwargs.items():
		args = BAKED_ARGS.get(k, {})
		if isinstance(v, dict):
			args.update(v)
		else:
			args.update({'path': str(v)})
		BAKED_ARGS[k] = args

import sh
# run command on foreground
sh_fg    = sh(_fg = True)
# piped command, waiting for piping
sh_piped = sh(_piped = True)
# redirect output
sh_out   = sh(_out = '>')
# piped and redirecting
sh_pout  = sh(_piped = True, _out = '>')
# run with Popen shell = True
sh_shell = sh(_shell = True)

sh.rm_rf = sh.rm.bake(r = True, f = True)
sh.ln_s  = sh.ln.bake(s = True)


