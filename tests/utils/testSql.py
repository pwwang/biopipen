import helpers, unittest, testly
from bioprocs.utils import sql

class TestUtilsSql(helpers.TestCase):

	def dataProvider_testDsnparse(self):
		yield 'sqlite:file=./x.db; encoding=utf8', dict(
			file     = './x.db',
			encoding = 'utf8',
			scheme   = 'sqlite',
			after    = ''
		)

	def testDsnparse(self, dsn, out):
		ret = sql.dsnparse(dsn)
		self.assertDictEqual(ret, out)

if __name__ == '__main__':
	testly.main(verbosity = 2)
