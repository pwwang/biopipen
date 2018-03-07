# PYPPL REPEAT START: readRecord
if 'readRecord' not in vars() or not callable(readRecord):
	class readRecord(object):

		def __repr__(self):
			return 'readRecord(%s)' % (', '.join([k + '=' + v for k, v in self.__dict__['__data'].items()]))

		def __init__(self, **kwargs):
			self.__dict__['__data'] = {}
			self.add(**kwargs)
		
		def keys(self):
			return self.__dict__['__data'].keys()

		def __getattr__(self, name):
			if name == '__data': return self.__dict__['__data']
			return self.__dict__['__data'][name]

		def __getitem__(self, name):
			return self.__getattr__(name)

		def __setitem__(self, name, val):
			self.__setattr__(name, val)
			
		def __setattr__(self, name, val):
			self.__dict__['__data'][name] = val
		
		def add(self, **kwargs):
			self.__dict__['__data'].update(kwargs)

		def haskey(self, key):
			return key in self.__dict__['__data']
# PYPPL REPEAT END: readRecord