import helpers, unittest
from os import path
from collections import OrderedDict
from bioprocs.utils.tsvio import TsvMeta, TsvRecord, TsvReader, TsvWriter, tsvOps, TsvioFiletypeException

class TestTsvMeta(helpers.TestCase):
	
	def dataProvider_testInit(self):
		yield TsvMeta(),
		
	def testInit(self, meta):
		self.assertIsInstance(meta, TsvMeta)
		self.assertIsInstance(meta, OrderedDict)
		self.assertDictEqual(meta, {})
		
	def dataProvider_testAdd(self):
		yield ['a'], {}, {'a': None}
		yield ['a', 'b'], {}, {'a': None, 'b': None}
		yield [], {'a': None}, {'a': None}
		yield [], {'a': 'b'}, {}, TypeError
		yield ['a', 'b', ['c', 'd']], {'e': int, 'f': str}, {
			'a': None,
			'b': None,
			'c': None,
			'd': None,
			'e': int,
			'f': str,
		}
		
	def testAdd(self, args, kwargs, data, exception = None):
		meta = TsvMeta()
		if exception:
			self.assertRaises(exception, meta.add, *args, **kwargs)
		else:
			meta.add(*args + kwargs.items())
			self.assertDictEqual(meta, data)	
			
	def dataProvider_testKeys(self):
		yield ['a', 'b', ['c', 'd']], {'e': int, 'f': str}, ['a', 'b', 'c', 'd', 'e', 'f']
	
	def testKeys(self, args, kwargs, keys):
		meta = TsvMeta(*args + kwargs.items())
		self.assertListEqual(meta.keys(), keys)
	
	def dataProvider_testRepr(self):
		yield TsvMeta(), 'TsvMeta()'
		yield TsvMeta('a', 'b', ('c', int)), 'TsvMeta(a, b, c=int)'
		yield TsvMeta('a', 'b', ('c', lambda x: int(x))), 'TsvMeta(a, b, c=<lambda>)'
	
	def testRepr(self, meta, reprstr):
		self.assertEqual(repr(meta), reprstr)
	
	def dataProvider_testGetAttr(self):
		meta = TsvMeta('a', 'b', ('c', int))
		yield meta, 'a', None
		yield meta, 'b', None
		yield meta, 'c', int
	
	def testGetAttr(self, meta, name, val):
		self.assertEqual(getattr(meta, name), val)
		
	def dataProvider_testSetAttr(self):
		meta = TsvMeta('a', 'b', ('c', int))
		yield meta, 'c', str
		yield meta, 'a', int
		yield meta, 'data', None
		
	def testSetAttr(self, meta, name, val):
		setattr(meta, name, val)
		self.assertEqual(getattr(meta, name), val)
	
	def dataProvider_testBorrow(self):
		meta1 = TsvMeta()
		meta2 = TsvMeta('a', 'b')
		yield meta1, meta2, {'a':None, 'b':None}
		meta1 = TsvMeta('a', 'b')
		meta2 = TsvMeta(('a', int), ('b', str))
		yield meta1, meta2, {'a':int, 'b':str}
	
	def testBorrow(self, meta1, meta2, data):
		meta1.update(meta2)
		self.assertDictEqual(meta1, data)
		
class TestTsvRecord(helpers.TestCase):
	
	def dataProvider_testInit(self):
		yield {'a':1, 'b':2, 'c':3}, {'a':1, 'b':2, 'c':3}
		
	def testInit(self, kwargs, data):
		record = TsvRecord(**kwargs)
		self.assertDictEqual(record, data)
		
	def dataProvider_testAdd(self):
		yield {'a':1, 'b':2, 'c':3}, {'a':1, 'b':2, 'c':3}
		
	def testAdd(self, kwargs, data):
		record = TsvRecord()
		record.update(**kwargs)
		self.assertDictEqual(record, data)
		
	def dataProvider_testKeys(self):
		yield {'a':1, 'b':2, 'c':3}, ['a', 'b', 'c']
		
	def testKeys(self, kwargs, keys):
		record = TsvRecord(**kwargs)
		self.assertItemEqual(record.keys(), keys)
		
	def dataProvider_testGetAttr(self):
		record = TsvRecord(a=1, b=2, c=3)
		yield record, 'a', 1
		yield record, 'b', 2
		yield record, 'c', 3
	
	def testGetAttr(self, record, name, val):
		self.assertEqual(getattr(record, name), val)
		self.assertEqual(record[name], val)
		
	def dataProvider_testSetAttr(self):
		record = TsvRecord(a=1, b=2, c=3)
		yield record, 'c', str
		yield record, 'a', int, False
		yield record, 'data', None, False
		
	def testSetAttr(self, record, name, val, setAttr = True):
		if setAttr:
			setattr(record, name, val)
		else:
			record[name] = val
		self.assertEqual(getattr(record, name), val)	
	
	def dataProvider_testSetAttr(self):
		record = TsvRecord(a=1, b=2, c=3)
		yield record, 'c', str
		yield record, 'a', int, False
		yield record, 'data', None, False
		
	def dataProvider_testHasKey(self):
		record = TsvRecord(a=1, b=2, c=3)
		yield record, 'a', True
		
	def testHasKey(self, record, key, has):
		self.assertEqual(key in record, has)
		
class TestTsvReaderBase(helpers.TestCase):

	def dataProvider_testInit(self, testdir, indir):
		yield path.join(indir, 'readerbase1.txt'), {}
		yield path.join(indir, 'readerbase1.txt'), {
			'delimit': ',',
			'comment': '#',
			'skip': 1
		}
	
	def testInit(self, infile, inopts):
		reader = TsvReader(infile, **inopts)
		self.assertEqual(type(reader).__name__, 'TsvReaderBase')
		self.assertIsInstance(reader.meta, TsvMeta)
		self.assertEqual(type(reader.file).__name__, 'file')
		self.assertEqual(reader.delimit, inopts['delimit'] if 'delimit' in inopts else '\t')
		self.assertEqual(reader.comment, inopts['comment'] if 'comment' in inopts else '#')
		self.assertEqual(reader.tell, 0)
	
	def dataProvider_testParse(self, testdir, indir):
		infile = path.join(indir, 'readerbase1.txt')
		reader = TsvReader(infile)
		reader.meta.add('A', 'B', ('C', int))
		yield reader, ['a', 'b', '1'], TsvRecord({'A':'a', 'B':'b', 'C':1})
	
	def testParse(self, reader, line, record):
		r = reader._parse(line)
		self.assertDictEqual(r, record)
		
	def dataProvider_testNext(self, testdir, indir):
		infile = path.join(indir, 'readerbase1.txt')
		yield TsvReader(infile), {}, StopIteration
		
		infile = path.join(indir, 'readerbase2.txt')
		reader = TsvReader(infile)
		reader.meta.add('A', 'B', ('C', int))
		yield reader, {'A':'a', 'B':'b', 'C':1}
		
	def testNext(self, reader, record, exception = None):
		if exception:
			self.assertRaises(exception, reader.next)
		else:
			r = reader.next()
			self.assertDictEqual(r, record)
	
	def dataProvider_testDump(self, testdir, indir):
		infile = path.join(indir, 'readerbase2.txt')
		reader = TsvReader(infile)
		reader.meta.add('A', 'B', ('C', int))
		yield reader, [{'A':'a', 'B':'b', 'C':1}, {'A':'a2', 'B':'b2', 'C':2}, {'A':'a3', 'B':'b3', 'C':3}]
	
	def testDump(self, reader, records):
		for i, r in enumerate(reader.dump()):
			self.assertDictEqual(r, records[i])
			
	def dataProvider_testRewind(self, testdir, indir):
		infile = path.join(indir, 'readerbase2.txt')		
		reader = TsvReader(infile)
		reader.meta.add('A', 'B', ('C', int))
		reader.next()
		yield reader,  {'A':'a', 'B':'b', 'C':1}
	
	def testRewind(self, reader, recaw):
		reader.rewind()
		self.assertDictEqual(reader.next(), recaw)
	
	def dataProvider_testAutoMeta(self, testdir, indir):
		yield path.join(indir, 'readerbase1.txt'), ['COL1']
		yield path.join(indir, 'readerbase2.txt'), ['COL1', 'COL2', 'COL3']
	
	def testAutoMeta(self, infile, meta):
		reader = TsvReader(infile)
		reader.autoMeta()
		self.assertListEqual(reader.meta.keys(), meta)
		
class TestTsvReaderBed(helpers.TestCase):
	
	def dataProvider_testInit(self, testdir, indir):
		infile = path.join(indir, 'readerbed1.txt')
		yield infile,
	
	def testInit(self, infile):
		reader = TsvReader(infile, ftype = 'bed')
		self.assertDictEqual(reader.meta, {
			'CHR'   : None,
			'START' : int,
			'END'   : int,
			'NAME'  : None,
			'SCORE' : float,
			'STRAND': None
		})
		self.assertEqual(reader.index, 1)
		
	def dataProvider_testParse(self, testdir, indir):
		infile = path.join(indir, 'readerbed1.txt')
		reader = TsvReader(infile, ftype = 'bed')
		yield reader, ['chr1', '1', '1000', 'seg1', '1000', '-']
		
	def testParse(self, reader, line):
		r = reader._parse(line)
		self.assertEqual(r.CHR, line[0])
		self.assertEqual(r.START, int(line[1]))
		self.assertEqual(r.END, int(line[2]))
		self.assertEqual(r.NAME, line[3])
		self.assertEqual(r.SCORE, float(line[4]))
		self.assertEqual(r.STRAND, line[5])
		
class TestTsvReaderBed12(helpers.TestCase):
	
	def dataProvider_testInit(self, testdir, indir):
		infile = path.join(indir, 'readerbed12.txt')
		yield infile,
	
	def testInit(self, infile):
		reader = TsvReader(infile, ftype = 'bed12')
		self.assertDictEqual(reader.meta, {
			'CHR'         : None,
			'START'       : int,
			'END'         : int,
			'NAME'        : None,
			'SCORE'       : float,
			'STRAND'      : None,
			'THICKSTART'  : int,
			'THICKEND'    : int,
			'ITEMRGB'     : None,
			'BLOCKCOUNT'  : int,
			'BLOCKSIZES'  : None,
			'BLOCKSTARTS' : None
		})
		
	def dataProvider_testParse(self, testdir, indir):
		infile = path.join(indir, 'readerbed12.txt')
		reader = TsvReader(infile, ftype = 'bed12')
		yield reader, ['chr1', '1', '1000', 'seg1', '1000', '-', '1', '100', '244', '8', '1,101', '100,201']
		
	def testParse(self, reader, line):
		r = reader._parse(line)
		self.assertEqual(r.CHR, line[0])
		self.assertEqual(r.START, int(line[1]))
		self.assertEqual(r.END, int(line[2]))
		self.assertEqual(r.NAME, line[3])
		self.assertEqual(r.SCORE, float(line[4]))
		self.assertEqual(r.STRAND, line[5])	
		self.assertEqual(r.THICKSTART, int(line[6]))
		self.assertEqual(r.THICKEND, int(line[7]))
		self.assertEqual(r.ITEMRGB, line[8])	
		self.assertEqual(r.BLOCKCOUNT, int(line[9]))
		self.assertEqual(r.BLOCKSIZES, line[10])	
		self.assertEqual(r.BLOCKSTARTS, line[11])	

class TestTsvReaderBedpe(helpers.TestCase):
	
	def dataProvider_testInit(self, testdir, indir):
		infile = path.join(indir, 'readerbedpe.txt')
		yield infile,
	
	def testInit(self, infile):
		reader = TsvReader(infile, ftype = 'bedpe')
		self.assertDictEqual(reader.meta, {
			'CHR1'    : None,
			'START1'  : int,
			'END1'    : int,
			'CHR2'    : None,
			'START2'  : int,
			'END2'    : int,
			'NAME'    : None,
			'SCORE'   : float,
			'STRAND1' : None,
			'STRAND2' : None
		})
		
	def dataProvider_testParse(self, testdir, indir):
		infile = path.join(indir, 'readerbedpe.txt')
		reader = TsvReader(infile, ftype = 'bedpe')
		yield reader, ['chr1', '1', '1000', 'chr1', '1000', '2000', 'interaction1', '100', '+', '-']
		
	def testParse(self, reader, line):
		r = reader._parse(line)
		self.assertEqual(r.CHR1, line[0])
		self.assertEqual(r.START1, int(line[1]))
		self.assertEqual(r.END1, int(line[2]))
		self.assertEqual(r.CHR2, line[3])
		self.assertEqual(r.START2, int(line[4]))
		self.assertEqual(r.END2, int(line[5]))
		self.assertEqual(r.NAME, line[6])
		self.assertEqual(r.SCORE, float(line[7]))
		self.assertEqual(r.STRAND1, line[8])	
		self.assertEqual(r.STRAND2, line[9])	
		
class TestTsvReaderBedx(helpers.TestCase):
	
	def dataProvider_testInit(self, testdir, indir):
		infile = path.join(indir, 'readerbedx.txt')
		yield infile,
	
	def testInit(self, infile):
		reader = TsvReader(infile, ftype = 'bedx')
		self.assertDictEqual(reader.meta, {
			'CHR'   : None,
			'START' : int,
			'END'   : int,
			'NAME'  : None,
			'SCORE' : float,
			'STRAND': None
		})
		
	def dataProvider_testParse(self, testdir, indir):
		infile = path.join(indir, 'readerbedx.txt')
		reader = TsvReader(infile, ftype = 'bedx', xcols = ['REF', 'ALT'])
		yield reader, ['chr1', '1', '1000', 'seg1', '1000', '-', 'A', 'C']
		
	def testParse(self, reader, line):
		r = reader._parse(line)
		self.assertEqual(r.CHR, line[0])
		self.assertEqual(r.START, int(line[1]))
		self.assertEqual(r.END, int(line[2]))
		self.assertEqual(r.NAME, line[3])
		self.assertEqual(r.SCORE, float(line[4]))
		self.assertEqual(r.STRAND, line[5])
		self.assertEqual(r.REF, line[6])
		self.assertEqual(r.ALT, line[7])
		
class TestTsvReaderHead(helpers.TestCase):
	
	def dataProvider_testInit(self, testdir, indir):
		infile = path.join(indir, 'readerhead.txt')
		yield infile, {'skip': 4, 'tmeta': {'C': int}}, {'ROWNAMES':None, 'B':None, 'C':int}, {'ROWNAMES':'a', 'B': 'b', 'C': 1}
		
	def testInit(self, infile, inopts, meta, record):
		reader = TsvReader(infile, ftype = 'head', **inopts)
		self.assertDictEqual(reader.meta, meta)
		r = reader.next()
		self.assertDictEqual(r, record)
	
class TestTsvWriterBase(helpers.TestCase):
	
	def dataProvider_testInit(self, testdir, indir, outdir):
		infile  = path.join(indir , 'writerbase.txt')
		outfile = path.join(outdir, 'writerbase.txt')
		yield infile, outfile
	
	def testInit(self, infile, outfile):
		writer = TsvWriter(infile)
		writer.meta.add('A', 'B', 'C')
		record = TsvRecord(A='a', B='b', C=1)
		writer.writeHead()
		writer.write(record)
		writer.close()
		self.assertFileEqual(infile, outfile)
		
class TestTsvWriterBed(helpers.TestCase):
	
	def dataProvider_testInit(self, testdir, indir, outdir):
		infile  = path.join(indir , 'writerbed.txt')
		outfile = path.join(outdir, 'writerbed.txt')
		yield infile, outfile
	
	def testInit(self, infile, outfile):
		writer = TsvWriter(infile, ftype = 'bed')
		record = TsvRecord(CHR='chr1', START=1, END=1000, NAME='seq1', SCORE=1, STRAND='-')
		writer.writeHead()
		writer.write(record)
		writer.close()
		self.assertFileEqual(infile, outfile)

class TestTsvOps(helpers.TestCase):
	
	def dataProvider_testTsvOps(self, testdir, indir, outdir):
		infile  = path.join(indir, 'tsvops.txt')
		outfile = path.join(testdir, 'tsvops.txt')
		exptfile  = path.join(outdir, 'tsvops1.txt')
		exptfile2  = path.join(outdir, 'tsvops2.txt')
		
		yield infile, outfile, {}, {}, None, exptfile, TsvioFiletypeException
		yield infile, outfile, {'ftype': 'bed'}, {}, None, infile
		yield infile, outfile, {'ftype': 'c1,c2,c3,c4,c5,c6'}, {}, None, infile
		yield infile, outfile, {'ftype': 'c1,c2,c3,c4,c5,c6'}, {'head': True, 'headTransform': lambda cols: [c if c!='c2' else c.upper() for c in cols]}, None, exptfile
		yield infile, outfile, {'ftype': 'c1,c2,c3,c4,c5,c6', 'skip': 8}, {
			'head': True, 
			'headPrefix': '##', 
			'ftype': ['c1', 'c2', 'c3'], 
			'headTransform': lambda cols: ['Chr' if c == 'c1' else 'Start' if c == 'c2' else 'End' for c in cols]
		}, lambda row: setattr(row, 'c3', int(row.c3) + 1) or row, exptfile2
	
	def testTsvOps(self, infile, outfile, inopts, outopts, transform, exptfile, exception = None):
		if exception:
			self.assertRaises(exception, tsvOps, infile, outfile, inopts, outopts, transform)
		else:
			tsvOps(infile, outfile, inopts, outopts, transform)
			self.assertFileEqual(outfile, exptfile)
	
if __name__ == '__main__':
	unittest.main(verbosity = 2)