import unittest
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.algorithm import pRWR, pAR

class testAlgorithm (unittest.TestCase):
	
	def testRWR (self):
		self.maxDiff = None

		pRWR.input = (getfile("test.rwr.w.txt"), getfile("test.rwr.e.txt"))
		PyPPL(config).start(pRWR).run()
		procOK(pRWR, "test.rwr.ret", self)

	def testAR(self):
		pAR1 = pAR.copy()
		pAR1.input         = (
			getfile('ar.d.txt.gz'), 
			getfile('ar.pt.txt.gz'), 
			getfile('ar.y.txt.gz')
		)
		pAR1.args.parallel = True
		pAR1.args.svdP     = 25
		PyPPL(config).start(pAR1).run()
		procOK(pAR1, "ar.w.txt", self)
	
	@unittest.skip('ADMM takes too long.')
	def testARADMM(self):
		pAR2 = pAR.copy()
		pAR2.input         = (
			getfile('ar.d.txt.gz'), 
			getfile('ar.pt.txt.gz'), 
			getfile('ar.y.txt.gz')
		)
		pAR2.args.parallel = False
		pAR2.args.method   = 'admm'
		PyPPL(config).start(pAR2).run()

		
if __name__ == '__main__':
	unittest.main()
		