"""
Reader and writer for tsv file.
"""
from bioprocs.utils import alwaysList
from collections import OrderedDict

__all__ = ['TsvMeta', 'TsvRecord', 'TsvReader', 'TsvWriter', 'tsvOps', 'TsvioFiletypeException']

class TsvMeta(OrderedDict):
	"""
	Tsv meta data
	"""
	def __init__(self, *args):
		"""
		arg could be a string, a tuple or a list:
		'a', ('b': int) or ['c', 'd']
		"""
		super(TsvMeta, self).__init__()
		self.add(*args)

	def __repr__(self):
		return 'TsvMeta(%s)' % ', '.join([k if not v else '%s=%s'%(k,v.__name__) for k,v in self.items()])

	def __getattr__(self, name):
		if not name.startswith('_OrderedDict'):
			return self[name]
		super(TsvMeta, self).__getattr__(name)
		
	def __setattr__(self, name, val):
		if not name.startswith('_OrderedDict'):
			self[name] = val
		else:
			super(TsvMeta, self).__setattr__(name, val)
	
	def add(self, *args):
		for arg in args:
			if isinstance(arg, list):
				for a in arg: self[a] = None
			elif isinstance(arg, tuple):
				if arg[1] and not callable(arg[1]):
					raise TypeError('Expect callable for meta value.')
				self[arg[0]] = arg[1]
			else:
				self[arg] = None	
		
class TsvRecord(dict):
		
	def __repr__(self):
		return 'TsvRecord(%s)' % (', '.join([k + '=' + repr(v) for k, v in self.items()]))
	
	def __getattr__(self, name):
		return self[name]

	def __setattr__(self, name, val):
		self[name] = val
		
class TsvReaderBase(object):
	def __init__(self, infile, delimit = '\t', comment = '#', skip = 0):
		openfunc = open
		if infile.endswith('.gz'):
			import gzip
			openfunc = gzip.open

		self.meta    = TsvMeta()
		self.file    = openfunc(infile)
		self.delimit = delimit
		self.comment = comment
		self.tell    = 0
		
		if skip > 0:
			for _ in range(skip):
				self.file.readline()
		self.tell = self.file.tell()
		
	def autoMeta(self, prefix = 'COL'):
		line = self.file.readline()
		while self.comment and line.startswith(self.comment):
			line = self.file.readline()
		self.rewind()
		line = line.rstrip('\n').split(self.delimit)
		cols = [prefix + str(i+1) for i in range(len(line))]
		self.meta.add(*cols)

	def _parse(self, line):
		record = TsvRecord()
		for i, key in enumerate(self.meta.keys()):
			try:
				record[key] = self.meta[key](line[i]) if self.meta[key] else line[i]
			except IndexError:
				record[key] = ''
		return record

	def next(self):
		line = self.file.readline()
		while self.comment and line.startswith(self.comment):
			line = self.file.readline()
		line = line.rstrip('\n')
		# empty lines not allowed
		if not line: raise StopIteration()
		return self._parse(line.split(self.delimit))

	def dump(self):
		return [r for r in self]

	def rewind(self):
		self.file.seek(self.tell)

	def __iter__(self):
		return self

	def __del__(self):
		self.close()
			
	def close(self):
		if self.file:
			self.file.close()
			
class TsvReaderBed(TsvReaderBase):
	META = [
		('CHR'   , None),
		('START' , int),
		('END'   , int),
		('NAME'  , None),
		('SCORE' , float),
		('STRAND', None)
	]

	def __init__(self, infile, skip = 0, comment = '#', delimit = '\t'):
		super(TsvReaderBed, self).__init__(infile, skip = skip, comment = comment, delimit = delimit)
		self.meta = TsvMeta(*TsvReaderBed.META)
		self.index = 1

	def _parse(self, line):
		r = super(TsvReaderBed, self)._parse(line)
		if not r.NAME: 
			r.NAME = 'BED' + str(self.index)
			self.index += 1
		if not r.SCORE: r.SCORE   = 0.0
		if not r.STRAND: r.STRAND = '+'
		return r
		
class TsvReaderBed12(TsvReaderBase):
	META = [
		('CHR'         , None),
		('START'       , int),
		('END'         , int),
		('NAME'        , None),
		('SCORE'       , float),
		('STRAND'      , None),
		('THICKSTART'  , int),
		('THICKEND'    , int),
		('ITEMRGB'     , None),
		('BLOCKCOUNT'  , int),
		('BLOCKSIZES'  , None),
		('BLOCKSTARTS' , None)
	]

	def __init__(self, infile, skip = 0, comment = '#', delimit = '\t'):
		super(TsvReaderBed12, self).__init__(infile, skip = skip, comment = comment, delimit = delimit)
		self.meta = TsvMeta(*TsvReaderBed12.META)
		
class TsvReaderBedpe(TsvReaderBase):
	META = [
		('CHR1'    , None),
		('START1'  , int),
		('END1'    , int),
		('CHR2'    , None),
		('START2'  , int),
		('END2'    , int),
		('NAME'    , None),
		('SCORE'   , float),
		('STRAND1' , None),
		('STRAND2' , None)
	]

	def __init__(self, infile, skip = 0, comment = '#', delimit = '\t'):
		super(TsvReaderBedpe, self).__init__(infile, skip = skip, comment = comment, delimit = delimit)
		self.meta = TsvMeta(*TsvReaderBedpe.META)

class TsvReaderBedx(TsvReaderBase):
	META = [
		('CHR'   , None),
		('START' , int),
		('END'   , int),
		('NAME'  , None),
		('SCORE' , float),
		('STRAND', None)
	]

	def __init__(self, infile, skip = 0, comment = '#', delimit = '\t', xcols = None, headprefix = ''):
		super(TsvReaderBedx, self).__init__(infile, skip = skip, comment = comment, delimit = delimit)
		self.meta = TsvMeta(*TsvReaderBedx.META)
		
		xmeta = OrderedDict()
		if not xcols:
			pass
		elif isinstance(xcols, list):
			for xcol in xcols:
				xmeta[xcol] = None
		elif isinstance(xcols, dict):
			for xcol, callback in xcols.items():
				if not callable(callback):
					raise TypeError('Expect callable for xcols values.')
				xmeta[xcol] = callback
		else:
			xmeta[xcols] = None
		self.meta.add(*xmeta.items())
		
class TsvReaderHead(TsvReaderBase):

	def __init__(self, infile, comment = '#', delimit = '\t', skip = 0, tmeta = None):
		super(TsvReaderHead, self).__init__(infile, skip = skip, comment = comment, delimit = delimit)
		self.meta = TsvMeta()
		
		header = self.file.readline().strip('#\t\n ').split(delimit)
		self.tell = self.file.tell()
		row1   = self.file.readline().strip().split(delimit)
		if len(row1) == len(header) + 1:
			header.insert(0, 'ROWNAMES')
		metatype = OrderedDict()
		for head in header:
			metatype[head] = None if not tmeta or not isinstance(tmeta, dict) or not head in tmeta or not callable(tmeta[head]) else tmeta[head]
		self.meta.add(*metatype.items())
		
		self.rewind()
		
class TsvWriterBase(object):
	def __init__(self, outfile, delimit = '\t'):
		openfunc = open
		if outfile.endswith('.gz'):
			import gzip
			openfunc = gzip.open
		
		self.delimit = delimit
		self.meta    = TsvMeta()
		self.file    = open(outfile, 'w')
	
	def writeHead(self, prefix = '', delimit = None, transform = None):
		delimit = delimit or self.delimit
		keys = self.meta.keys()
		if callable(transform): 
			keys = transform(keys)
		elif isinstance(transform, dict):
			keys = [key if not key in transform or not callable(transform[key]) else transform[key](key)]
		self.file.write(prefix + delimit.join(keys) + '\n')

	def write(self, record, delimit = None):
		delimit = delimit or self.delimit
		outs = []
		for key in self.meta.keys():
			outs.append(str(record[key]))
		self.file.write(delimit.join(outs) + '\n')

	def __del__(self):
		self.close()
			
	def close(self):
		if self.file:
			self.file.close()
		
class TsvWriterBed(TsvWriterBase):
		
	def __init__(self, outfile):
		super(TsvWriterBed, self).__init__(outfile)
		self.meta = TsvMeta(*TsvReaderBed.META)
		
def TsvReader(infile, **inopts):
	if not 'ftype' in inopts or not inopts['ftype']:
		if 'ftype' in inopts: del inopts['ftype']
		return TsvReaderBase(infile, **inopts)
	else:
		ftype = inopts['ftype']
		del inopts['ftype']
		ftype = ftype[0].upper() + ftype[1:].lower()
		klass = globals()['TsvReader' + ftype]
		return klass(infile, **inopts)
		

def TsvWriter(outfile, **outopts):
	if not 'ftype' in outopts:
		return TsvWriterBase(outfile, delimit = outopts['delimit'] if 'delimit' in outopts else '\t')
	else:
		ftype = outopts['ftype']
		del outopts['ftype']
		ftype = ftype[0].upper() + ftype[1:].lower()
		klass = globals()['TsvWriter' + ftype]
		return klass(outfile, **outopts)
	

class TsvioFiletypeException(Exception):
	pass

def tsvOps(infile, outfile, inopts, outopts, transform = None):
	if 'ftype' not in inopts:
		raise TsvioFiletypeException('No file type specified.')
	try:
		reader = TsvReader(infile, **inopts)
	except (KeyError, TypeError):
		ftype = inopts['ftype']
		del inopts['ftype']
		reader = TsvReader(infile, **inopts)
		reader.meta.add(*alwaysList(ftype))
	
	try:
		writer = TsvWriter(outfile, **outopts)
		writer.meta.update(reader.meta)
	except (KeyError, TypeError, AttributeError):
		ftype = outopts['ftype']
		del outopts['ftype']
		writer = TsvWriter(outfile, **outopts)
		writer.meta.add(*alwaysList(ftype))
	
	if 'head' in outopts and outopts['head']:
		del outopts['head']
		if 'headPrefix' in outopts:
			headPrefix    = outopts['headPrefix']
			del outopts['headPrefix']
		else:
			headPrefix    = ''
		
		if 'headDelimit' in outopts:
			headDelimit   = outopts['headDelimit']
			del outopts['headDelimit']
		else:
			headDelimit   = None
			
		if 'headTransform' in outopts:
			headTransform = outopts['headTransform']
			del outopts['headTransform']
		else:
			headTransform = None
		writer.writeHead(headPrefix, headDelimit, headTransform)
	
	for r in reader:
		if callable(transform):
			r2 = transform(r)
			if not r2: continue
			writer.write(r2)
		else:
			writer.write(r)