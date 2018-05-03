import helpers, unittest, testly
from os import path
from bioprocs.utils.cache import Cache
from medoo import Medoo, MedooSqlite, Raw, Function

class TestCache(helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestCache')

	def dataProvider_testInit(self):
		dbfile = path.join(self.testdir, 'cache.db')
		yield dbfile, 'table1', {'id': 'int primary key', 'content': 'text'}, 'id'

	def testInit(self, dbfile, table, schema, pkey):
		c = Cache(dbfile, table, schema, pkey)
		self.assertIsInstance(c, Cache)
		self.assertIsInstance(c.medoo, MedooSqlite)
		self.assertEqual(c.table, table)
		self.assertTrue(c.medoo.tableExists(table))

	def dataProvider_testDummy(self):
		# 0
		yield 'plain', 'query', ('col', 'data'), 'col', 'data'
		yield 'plain', 'find', [{'col': 'qitem'}], 'col', 'qitem', [{'col': 'qitem'}]
		yield 'plain', 'find', [], 'col', 'qitem', [{'col': 'Qitem'}]
		yield 'plain', 'insert', ('col', 'data'), 'col', 'data'
		yield 'plain', 'update', ('col', 'data'), 'col', 'data', ''
		# 5
		yield 'iplain', 'query', (Function.upper('col[ IN ]'), ['DATA1', 'DATA2']), 'col', ['data1', 'data2']
		yield 'iplain', 'find', [{'col': 'qitem'}], 'col', 'qitem', [{'col': 'qitem'}]
		yield 'iplain', 'find', [{'col': 'Qitem'}], 'col', 'qitem', [{'col': 'Qitem'}]
		yield 'iplain', 'insert', ('col', 'data'), 'col', 'data'
		yield 'iplain', 'update', ('col', 'data'), 'col', 'data', ''
		# 10
		yield 'array', 'query', ('col[~]', ['% // data', '% // data // %']), 'col', ['data']
		yield 'array', 'find', [{'col': ['qitem', 'item2']}], 'col', 'qitem', [{'col': ['qitem', 'item2']}]
		yield 'array', 'find', [], 'col', 'qitem', [{'col': ' // Qitem // item2'}]
		yield 'array', 'insert', ('col', ' // data1 // data2'), 'col', ['data1', 'data2']
		yield 'array', 'update', ('col', '"col"||\' // data1 // data2\''), 'col', ['data1', 'data2'], ''
		# 15
		yield 'iarray', 'query', ('col[~~]', ['% // data', '% // data // %']), 'col', ['data']
		yield 'iarray', 'find', [{'col': ['qitem', 'item2']}], 'col', 'qitem', [{'col': ['qitem', 'item2']}]
		yield 'iarray', 'find', [{'col': ['Qitem', 'item2']}], 'col', 'qitem', [{'col': ['Qitem', 'item2']}]
		yield 'iarray', 'insert', ('col', ' // data1 // data2'), 'col', ['data1', 'data2']
		yield 'iarray', 'update', ('col', '"col"||\' // data1 // data2\''), 'col', ['data1', 'data2'], ''

	def testDummy(self, key, func, output, *args, **kwargs):
		if not isinstance(output, tuple):
			self.assertEqual(Cache.DUMMY[key][func](*args, **kwargs), output)
		else:
			self.assertListEqual(
				[str(o) for o in Cache.DUMMY[key][func](*args, **kwargs)],
				[str(o) for o in output]
			)

	def dataProvider_testQueryWhere(self):
		dummies = {'a':'plain', 'b':'iarray', 'c':{'query': lambda col, data: (col, col)}}
		yield ['a', 'b'], [1,2,3,4], dummies, {'OR#a:b': {'a': [1,2,3,4], 'b[~~]':['% // 1', '% // 1 // %', '% // 2', '% // 2 // %', '% // 3', '% // 3 // %', '% // 4', '% // 4 // %']}}
		yield 'c', [1,2,3,4], dummies, {'c': 'c'}

	def testQueryWhere(self, keys, data, dummies, output):
		self.maxDiff = None
		self.assertDictEqual(
			Cache._queryWhere(keys, data, dummies),
			output
		)

	def dataProvider_testFind(self):
		dummies1 = {'a':'plain', 'b':'iarray', 'c':{'find': lambda col, qitem, results: []}}
		dummies2 = {'a':'plain', 'b':'array', 'c':{'find': lambda col, qitem, results: []}}
		yield {'a, b': 'qitem'}, [{'a': '1', 'b': ['Qitem']}], dummies1, {'a': '1', 'b': ['Qitem']}
		yield {'a, b': 'qitem'}, [{'a': '1', 'b': ['Qitem']}], dummies2, None
		yield {'c': 'qitem'}, [{'c': ['123']}], dummies2, None

	def testFind(self, qitems, results, dummies, output):
		self.assertEqual(
			Cache._find(qitems, results, dummies),
			output
		)

	def dataProvider_testCheckData(self):
		yield {'a': 1, 'b': 2}, {'a': [1], 'b': [2]}
		yield {'a': [1,2], 'b': 3}, {'a': [1,2], 'b': [3,3]}
		yield {'a': [1,2], 'b': [1,2,3]}, {}, ValueError

	def testCheckData(self, indata, outdata, exception = None):
		if exception:
			self.assertRaises(exception, Cache._checkData, indata)
		else:
			l = Cache._checkData(indata)
			self.assertEqual(l, len(outdata.values()[0]))
			self.assertDictEqual(indata, outdata)

	'''
	def dataProvider_testItem2where(self):
		yield 'id', [1,2,3,4], {}, {'id': [1,2,3,4]}
		yield 'id', [1,2,3,4], {'id': 'id[!]'}, {'id[!]': [1,2,3,4]}
		yield 'id', [1,2,3,4], {'id': ('id[>]', lambda x: x+1)}, {
			'OR # id': {'id[>] # 0': 2, 'id[>] # 1': 3, 'id[>] # 2': 4, 'id[>] # 3': 5}
		}

	def testItem2where(self, key, val, where, ret):
		r = Cache._item2where(key, val, where)
		self.assertDictEqual(r, ret)
	'''
	def dataProvider_testQueryN(self):
		dbfile = ':memory:'
		c = Cache(dbfile, 'query', {
			'id': 'int',
			'c1': 'text',
			'c2': 'text'
		}, 'id')
		c.medoo.insert('query', {
			'id': 1, 'c1': 'G1', 'c2': ' // A1 // ALIAS1 // ANAME1'
		}, {
			'id': 2, 'c1': 'G2', 'c2': ' // A2 // ALIAS2 // ANAME2'
		}, {
			'id': 3, 'c1': 'G3', 'c2': ' // A3 // ALIAS3 // ANAME3'
		})
		c.medoo.insert('query', {
			'id': 4, 'c2': ' // A4 // ALIAS4 // ANAME4'
		})
		columns1 = ['c1']
		data1 = {
			'c2': ['A1', 'ALiAS2', 'ANAME4', 'ANAMe5', 'aname4']
		}
		dummies = {
			'c2': 'iarray'
		}
		datafound1 = {
			0: {'id': 1, 'c1': 'G1', 'c2': ['A1', 'ALIAS1', 'ANAME1']},
			1: {'id': 2, 'c1': 'G2', 'c2': ['A2', 'ALIAS2', 'ANAME2']}
		}
		datarest1 = {
			'c2': ['ANAME4', 'ANAMe5', 'aname4']
		}
		yield testly.Data(
			c, columns1, data1, dummies, datafound1, datarest1
		)


	def testQueryN(self, cache, columns, data, dummies, datafound, datarest):
		self.maxDiff = None
		retfound, retrest = cache._queryN(columns, data, dummies)
		self.assertDictEqual(retfound, datafound)
		self.assertItemsEqual(retrest.keys(), datarest.keys())
		for key in retrest.keys():
			self.assertItemsEqual(retrest[key], datarest[key])

	def dataProvider_testQuery(self):
		dbfile = ':memory:'
		c = Cache(dbfile, 'query', {
			'id': 'int',
			'c1': 'text',
			'c2': 'text'
		}, 'id')
		c.medoo.insert('query', {
			'id': 1, 'c1': 'G1', 'c2': ' // A1 // ALIAS1 // ANAME1'
		}, {
			'id': 2, 'c1': 'G2', 'c2': ' // A2 // ALIAS2 // ANAME2'
		}, {
			'id': 3, 'c1': 'G3', 'c2': ' // A3 // ALIAS3 // ANAME3'
		})
		columns1 = ['c1']
		data1 = {
			'c2': ['A1', 'ALIAS2', 'ANAME4', 'ANAME5']
		}
		dummies = {
			'c2': 'array'
		}
		datafound1 = {
			0: {'id': 1, 'c1': 'G1', 'c2': ['A1', 'ALIAS1', 'ANAME1']},
			1: {'id': 2, 'c1': 'G2', 'c2': ['A2', 'ALIAS2', 'ANAME2']}
		}
		datarest1 = {
			'c2': ['ANAME4', 'ANAME5']
		}
		yield c, columns1, data1, dummies, 1, datafound1, datarest1
		yield c, columns1, data1, dummies, 2, datafound1, datarest1
		yield c, columns1, data1, dummies, 3, datafound1, datarest1
		yield c, columns1, data1, dummies, 4, datafound1, datarest1
		yield c, columns1, data1, dummies, 5, datafound1, datarest1

	def testQuery(self, cache, columns, data, dummies, chunk, datafound, datarest):
		self.maxDiff = None
		retfound, retrest = cache.query(columns, data, dummies, chunk)
		self.assertDictEqual(retfound, datafound)
		self.assertItemsEqual(retrest.keys(), datarest.keys())
		for key in retrest.keys():
			self.assertItemsEqual(retrest[key], datarest[key])

	def dataProvider_testSave(self):
		dbfile = ':memory:'
		c = Cache(dbfile, 'query', {
			'id': 'int',
			'c1': 'text',
			'c2': 'text'
		}, 'id')
		c.medoo.insert('query', {
			'id': 4, 'c1': 'G4', 'c2': ' // A4'
		}, {
			'id': 5, 'c1': 'G5', 'c2': ''
		})
		yield c, {
			'id': [4, 5, 6],
			'c2': ['ANAME4', 'ANAME5', 'A6']
		}, {
			'c2': 'array'
		}, [
			{'id': 4, 'c1': 'G4', 'c2': ' // A4 // ANAME4'},
			{'id': 5, 'c1': 'G5', 'c2': ' // ANAME5'},
			{'id': 6, 'c1': None, 'c2': ' // A6'},
		]

	def testSave(self, cache, data, dummies, datall):
		self.maxDiff = None
		cache.save(data, dummies)
		rs = cache.medoo.select(cache.table, '*')
		self.assertListEqual(rs.fetchall(), datall)

if __name__ == '__main__':
	testly.main(verbosity = 2)
