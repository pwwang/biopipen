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
	picard = dict(_sep = '=', _prefix = '')
))

def _update_args(*args, **kwargs):
	for arg in args:
		assert isinstance(arg, dict)
		_update_args(**arg)
	for k, v in kwargs.items():
		args = BAKED_ARGS.get(k, {})
		if isinstance(v, dict):
			args.update(v)
		else:
			args.update({'path': str(v)})
		BAKED_ARGS[k] = args

import sh
shell = sh




