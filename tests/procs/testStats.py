import unittest
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.stats import pMetaPval, pMetaPval1, pSurvival
from bioprocs.common import pFiles2Dir
from bioprocs.matrix import pRbind

unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)
class testStats (unittest.TestCase):
	def test1MetaPval(self):
		pMetaPval.input        = [getfile()]
		pMetaPval.args.pattern = "testMetaPval*.txt"
		PyPPL(config).start(pMetaPval).run()
		procOK(pMetaPval, 'metapval.txt', self)

	def test2MetaPval1(self):
		pRbind.input = [[getfile('testMetaPval%s.txt' % i) for i in range(10)]]
		pRbind.args.rnames = [False] * 10
		pRbind.args.header = True

		pMetaPval1.depends = pRbind
		
		PyPPL(config).start(pRbind).run()
		procOK(pMetaPval1, 'metapval2.txt', self)

	def testSurvival(self):
		pSurvival1 = pSurvival.copy()
		pSurvival1.input                = [getfile('survival.txt')]
		pSurvival1.args.rnames          = False
		pSurvival1.args.nthread         = 3
		pSurvival1.args.gridParams.nrow = 2
		pSurvival1.args.gridParams.ncol = 2
		PyPPL(config).start(pSurvival1).run()
		procOK(pSurvival1, 'survival-combine', self)


	def testSurvivalNoCombine(self):
		pSurvival2              = pSurvival.copy()
		pSurvival2.input        = [getfile('survival.txt')]
		pSurvival2.args.rnames  = False
		pSurvival2.args.nthread = 1
		pSurvival2.args.combine = False
		PyPPL(config).start(pSurvival2).run()
		procOK(pSurvival2, 'survival-individual', self)


if __name__ == '__main__':
	unittest.main(failfast = True)