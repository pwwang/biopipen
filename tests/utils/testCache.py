import helpers, unittest
from os import path
from bioprocs.utils.cache import Cache
from medoo import Medoo, MedooSqlite, Raw

class TestCache(helpers.TestCase):
	
	def dataProvider_testInit(self, testdir):
		dbfile = path.join(testdir, 'cache.db')
		yield dbfile, 'table1', {'id': 'int primary key', 'content': 'text'}, 'id'

	def testInit(self, dbfile, table, schema, pkey):
		c = Cache(dbfile, table, schema, pkey)
		self.assertIsInstance(c, Cache)
		self.assertIsInstance(c.medoo, MedooSqlite)
		self.assertEqual(c.table, table)
		self.assertTrue(c.medoo.tableExists(table))
	
	def dataProvider_testItem2where(self):
		yield 'id', [1,2,3,4], {}, {'id': [1,2,3,4]}
		yield 'id', [1,2,3,4], {'id': 'id[!]'}, {'id[!]': [1,2,3,4]}
		yield 'id', [1,2,3,4], {'id': ('id[>]', lambda x: x+1)}, {
			'OR # id': {'id[>] # 0': 2, 'id[>] # 1': 3, 'id[>] # 2': 4, 'id[>] # 3': 5}
		}		
	
	def testItem2where(self, key, val, where, ret):
		r = Cache._item2where(key, val, where)
		self.assertDictEqual(r, ret)
		
	def dataProvider_testQueryN(self, testdir):
		dbfile = ':memory:'
		c = Cache(dbfile, 'query', {
			'id': 'int',
			'c1': 'text',
			'c2': 'text'
		}, 'id')
		c.medoo.insert('query', {
			'id': 1, 'c1': 'G1', 'c2': 'A1|ALIAS1|ANAME1'
		}, {
			'id': 2, 'c1': 'G2', 'c2': 'A2|ALIAS2|ANAME2'
		}, {
			'id': 3, 'c1': 'G3', 'c2': 'A3|ALIAS3|ANAME3'
		})
		columns1 = ['c1']
		data1 = {
			'c2': ['A1', 'ALIAS2', 'ANAME4', 'ANAME5']
		}
		where1 = {
			'c2': ('c2[~]', lambda x: ['{}|%'.format(x), '%|{}'.format(x), '%|{}|%'.format(x)])
		}
		found1 = {
			'c2': lambda v, data: data.startswith('%s|' % v) or data.endswith('|%s' % v) or '|%s|' % v in data
		}
		datafound1 = [
			{'id': 1, 'c1': 'G1', 'c2': 'A1|ALIAS1|ANAME1'},
			{'id': 2, 'c1': 'G2', 'c2': 'A2|ALIAS2|ANAME2'}
		]
		datarest1 = {
			'c2': ['ANAME4', 'ANAME5']
		}
		yield c, columns1, data1, where1, found1, datafound1, datarest1
		
		
	def testQueryN(self, cache, columns, data, where, found, datafound, datarest):
		retfound, retrest = cache.queryN(columns, data, where, found)
		self.assertListEqual(retfound, datafound)
		self.assertItemsEqual(retrest.keys(), datarest.keys())
		for key in retrest.keys():
			self.assertItemsEqual(retrest[key], datarest[key])
		
	def dataProvider_testQuery(self, testdir):
		dbfile = ':memory:'
		c = Cache(dbfile, 'query', {
			'id': 'int',
			'c1': 'text',
			'c2': 'text'
		}, 'id')
		c.medoo.insert('query', {
			'id': 1, 'c1': 'G1', 'c2': 'A1|ALIAS1|ANAME1'
		}, {
			'id': 2, 'c1': 'G2', 'c2': 'A2|ALIAS2|ANAME2'
		}, {
			'id': 3, 'c1': 'G3', 'c2': 'A3|ALIAS3|ANAME3'
		})
		columns1 = ['c1']
		data1 = {
			'c2': ['A1', 'ALIAS2', 'ANAME4', 'ANAME5']
		}
		where1 = {
			'c2': ('c2[~]', lambda x: ['{}|%'.format(x), '%|{}'.format(x), '%|{}|%'.format(x)])
		}
		found1 = {
			'c2': lambda v, data: data.startswith('%s|' % v) or data.endswith('|%s' % v) or '|%s|' % v in data
		}
		datafound1 = [
			{'id': 1, 'c1': 'G1', 'c2': 'A1|ALIAS1|ANAME1'},
			{'id': 2, 'c1': 'G2', 'c2': 'A2|ALIAS2|ANAME2'}
		]
		datarest1 = {
			'c2': ['ANAME4', 'ANAME5']
		}
		yield c, columns1, data1, where1, found1, 1, datafound1, datarest1
		yield c, columns1, data1, where1, found1, 2, datafound1, datarest1
		
	def testQuery(self, cache, columns, data, where, found, chunk, datafound, datarest):
		retfound, retrest = cache.query(columns, data, where, found, chunk)
		self.assertListEqual(retfound, datafound)
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
			'id': 4, 'c1': 'G4', 'c2': 'A4'
		}, {
			'id': 5, 'c1': 'G5', 'c2': ''
		})
		yield c, {
			'id': [4, 5, 6],
			'c2': ['ANAME4', 'ANAME5', 'A6']
		}, {
			'c2': ((lambda v: Raw('"c2" || \'|%s\'' % v)), lambda v:v)
		}, [
			{'id': 4, 'c1': 'G4', 'c2': 'A4|ANAME4'},
			{'id': 5, 'c1': 'G5', 'c2': '|ANAME5'},
			{'id': 6, 'c1': None, 'c2': 'A6'},
		]
			
	def testSave(self, cache, data, factory, datall):
		cache.save(data, factory)
		rs = cache.medoo.select(cache.table, None, '*')
		self.assertListEqual(rs.all(), datall)
	
if __name__ == '__main__':
	unittest.main(verbosity = 2)