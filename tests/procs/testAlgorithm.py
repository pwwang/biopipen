import unittest
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.algorithm import pRWR

class testAlgorithm (unittest.TestCase):
	
	def testRWR (self):
		self.maxDiff = None

		pRWR.input = (getfile("test.rwr.w.txt"), getfile("test.rwr.e.txt"))
		PyPPL(config).start(pRWR).run()
		procOK(pRWR, "test.rwr.ret", self)
		
if __name__ == '__main__':
	unittest.main()
		