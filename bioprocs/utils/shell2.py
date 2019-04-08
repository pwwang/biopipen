# requires https://github.com/pwwang/sh
# install: pip install git+https://github.com/pwwang/sh
from pyppl import Box
from sh import BAKED_ARGS
# defaults
BAKED_ARGS.update(dict(
	sort = dict(_sep = '', _duplistkey = True)
))

def updatePaths(*args, **kwargs):
	for arg in args:
		assert isinstance(arg, dict)
		updatePaths(**arg)
	for k, v in kwargs.items():
		args = BAKED_ARGS.get(k, {})
		args.update({'path': str(v)})
		BAKED_ARGS[k] = args

update_paths = updatePaths

import sh
shell = sh




