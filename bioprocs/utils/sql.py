from os import path

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
	
def schemaparse(schemafile, type = 'schema', scheme = 'sqlite', delimit = '\t'):
	name, fields = None, []
	name = path.basename(schemafile).split('.')[0]
	if type == 'schema':
		with open(schemafile) as f:
			for line in f:
				line = line.strip(' \n')
				if not line: continue
				if line.startswith('#'):
					line = line.strip()
					if line.startswith('Name:'):
						name = line[5:].strip()
					else:
						continue
				parts = line.split(delimit)
				if len(parts) == 2: parts.append('')
				fields.append(parts)
	else: # type == 'data'
		if scheme not in ['sqlite', 'sqlite3']:
			raise NotImplementedError('Only support sqlite3 currently.')
		with open(schemafile) as f:
			header = f.readline().strip('# \n').split(delimit)
		for head in header:
			fields.append([head, 'TEXT', ''])
	return name, fields