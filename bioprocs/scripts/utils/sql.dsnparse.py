if 'dsnparse' not in vars() or not callable(dsnparse):
	def dsnparse(dsn):
		ret = lambda: None
		scheme, args = dsn.split(':', 1)
		setattr(ret, 'scheme', scheme)
		args = args.split(';')
		for arg in args:
			arg = arg.strip()
			if not arg: continue
			argname, argval = arg.split('=')
			setattr(ret, argname, argval)
		if not hasattr(ret, 'after'):
			setattr(ret, 'after', '')
		return ret