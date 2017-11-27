import unittest
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.stats import pMetaPval, pMetaPval1
from bioprocs.common import pFiles2Dir
from bioprocs.matrix import pRbind

unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)
class testStats (unittest.TestCase):
	def test1MetaPval(self):
		pMetaPval.input = [getfile()]
		PyPPL().start(pMetaPval).run()
		procOK(pMetaPval, 'metapval.txt', self)

	def test2MetaPval1(self):
		pRbind.input = [[getfile('testMetaPval%s.txt' % i) for i in range(10)]]
		pRbind.args.rnames = [False] * 10
		pRbind.args.header = True

		pMetaPval1.depends = pRbind
		
		PyPPL().start(pRbind).run()
		procOK(pMetaPval1, 'metapval2.txt', self)

if __name__ == '__main__':
	unittest.main(failfast = True)