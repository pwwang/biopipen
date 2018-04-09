from medoo import Medoo

class Cache(object):
	
	def __init__(self, dbfile, table, schema, pkey):
		# if dbfile exists, check table exists, check cols
		# else create the table
		
		self.medoo = Medoo(
			database_type = 'sqlite', 
			database_file = dbfile, 
			check_same_thread = True,
			logging = True
		)
		self.table = table
		self.pkey  = pkey
		self.medoo.create(table, schema)
	
	@staticmethod
	def _item2where(key, val, where):
		conditions = {}
		#{'id': ([1,2,3,4], )}
		if not key in where:
			# 'id' IN (1,2,3,4)
			conditions[key] = val
		else:
			how = where[key]
			if not isinstance(how, tuple):
				how = (how, )
			if len(how) == 1:
				conditions[how[0]] = val
			else:
				conditions  = {'OR # %s' % key: {}}
				values = [
					how[1](v) for v in val
				]
				for i, v in enumerate(values):
					conditions['OR # %s' % key][how[0] + ' # ' + str(i)] = v
		return conditions
	
	"""
	data: {
		'id': [1,2,3,4]
		'alias': ['AT', 'AF', 'ATF', 'ATF1']
	}
	where: {
		'alias': ('alias[REGEXP]', lambda x: '^'+x+'$')
	}
	found: {
		'id': lambda v, data: v in data
		'alias': lambda v, data: any([v in d for d in data])
	}
	"""
	def queryN(self, columns, data, where, found = None):
		conditions = []
		for key, val in data.items():
			conditions.append(Cache._item2where(key, val, where))
		if len(conditions) == 1:
			where = conditions[0]
		else:
			where = {'AND': {}}
			for c in conditions:
				where['AND'].update(c)
		if not isinstance(columns, list):
			columns = [c.strip() for c in columns.split(',')]
		columns += data.keys() + [self.pkey]
		columns = list(set(columns))
		rs = self.medoo.select(self.table, columns, where)
		datall = rs.fetchall() if rs else []
		datafound = {}
		for key, val in data.items():
			foundfunc = found[key] if key in found else lambda v, d: v == d
			datafound[key] = [v for v in val if any([foundfunc(v, d[key]) for d in datall])]
		rest = {}
		for key, val in data.items():
			rest[key] = list(set(data[key]) - set(datafound[key]))
		
		return datall, rest
		
	def query(self, columns, data, where, found = None, chunk = 1000):
		retall, retrest = [], {}
		for i in xrange(0, len(data.values()[0]), chunk):
			datai = {key:val[i:i+chunk] for key, val in data.items()}
			datall, rest = self.queryN(columns, datai, where, found)
			retall += datall
			retrest = {key: (val + (retrest[key] if key in retrest else [])) for key, val in rest.items()}
		return retall, retrest
	"""
	data: {
			'id': [4, 5, 6],
			'c2': ['ANAME4', 'ANAME5', 'A6']
	}
	factory: {
		#         update                               insert
		'c2': ((lambda v: Raw('"c2" || \'|%s\'' % v)), lambda v:v)
	}
	"""
	def save(self, data, factory = None):
		factory = factory or {}
		pkdata  = data[self.pkey]
		rs      = self.medoo.select(self.table, self.pkey, {self.pkey: pkdata})
		rsall   = rs.fetchall() if rs else []
		dtidxs  = []
		funcup  = {}
		funcin  = {}
		for key, val in factory.items():
			if isinstance(val, (tuple, list)):
				if len(val) == 1:
					funcup[key] = val[0]
				elif len(val) > 1:
					funcup[key] = val[0]
					funcin[key] = val[1]
			else:
				funcup[key] = val
		for r in rsall:
			pkval  = r[self.pkey]
			try: 
				dtidx  = data[self.pkey].index(pkval)
			except ValueError:
				dtidx  = data[self.pkey].index(str(pkval))
			dtidxs.append(dtidx)
			updata = {}
			for key, val in data.items():
				if key == self.pkey: continue
				updata[key] = funcup[key](val[dtidx]) if key in funcup else val[dtidx]
				#del data[key][dtidx]
			self.medoo.update(self.table, updata, {self.pkey: pkval}, commit = False)
		
		self.medoo.commit()
		
		indata = []
		keys = data.keys()
		for i in range(len(data[keys[0]])):
			if i in dtidxs: continue
			indata.append({
				key:(funcin[key](data[key][i]) if key in funcin else data[key][i]) \
				for key in keys
			})
		if indata:
			self.medoo.insert(self.table, *indata)