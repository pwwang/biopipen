if 'readMeta' not in vars() or not callable(readMeta):
	from collections import OrderedDict
	class readMeta(object):
		def __init__(self, *args, **kwargs):
			"""
			readMeta((key1, val1), (key2, val2)) or
			readMeta(key1, key2)
			"""
			self.__dict__['__data'] = OrderedDict()
			self.add(*args, **kwargs)
		
		def keys(self):
			return self.__dict__['__data'].keys()

		def __repr__(self):
			return 'readMeta(%s)' % ', '.join(self.keys())

		def __getattr__(self, name):
			if name == '__data': return self.__dict__['__data']
			return self.__dict__['__data'][name]
			
		def __setattr__(self, name, val):
			self.__dict__['__data'][name] = val
		
		def borrow(self, x):
			self.__dict__['__data'].update(x.__dict__['__data'])
		
		def add(self, *args, **kwargs):
			for item in args:
				if isinstance(item, list) or isinstance(item, tuple):
					if len(item) == 0:
						raise ValueError('Expect at least 1 element.')
					elif len(item) == 1:
						key = item[0]
						val = ''
					else:
						key, val = item
				else:
					key = str(item)
					val = ''
				self.__dict__['__data'][key] = val
			if len(kwargs) > 1:
				raise ValueError('kwargs can only accept 1 key, as order will lost for multiple keys.')
			elif len(kwargs) == 1:
				key, val = kwargs.items()[0]
				self.__dict__['__data'][key] = val
		
		def remove(self, key):
			del self.__dict__['__data'][key]

