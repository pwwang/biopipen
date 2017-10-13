if 'params2CmdArgs' not in vars() or not callable(params2CmdArgs):
	# dash  = 'auto', '-', '--', ''
	# equal = 'auto', '=', ' '
	def params2CmdArgs(params, dash = 'auto', equal = 'auto', noq = None):
		try:  # py3
			from shlex import quote
		except ImportError:  # py2
			from pipes import quote
		ret = []
		for key, val in params.items():
			key = key.strip()
			if isinstance(val, bool) and not val: continue
			item = '--' if (dash == 'auto' and len(key) > 1) else '-' if dash == 'auto' else dash
			item += key
			if not isinstance(val, bool):
				item += '=' if (equal == 'auto' and len(key)>1) else ' ' if equal == 'auto' else equal
				if not noq or key not in list(noq):
					item += quote(str(val))
				else:
					item += str(val)
			ret.append(item)
		return ' '.join(ret)
