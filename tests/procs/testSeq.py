import unittest
from os import path
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs import params
from bioprocs.seq import pPromoters, pConsvPerm, pConsv
from bioprocs.bed import pBedRandom
from bioprocs.common import pStr2File

unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)

class TestSeq (unittest.TestCase):

	def testpPromoters(self):
		pStr2File.input = ['TP53']
		pPromoters.depends = pStr2File
		PyPPL(config).start(pStr2File).run()
		procOK(pPromoters, 'promoters.bed', self)

	@unittest.skipIf(not path.exists(params.consvdir.value), 'BigWig dir not exists.')
	def testpConsv(self):
		pBedRandom.input  = [(50, 100)]
		pConsvPerm.input  = [0]
		pConsv.depends    = pBedRandom, pConsvPerm
		pConsv2           = pConsv.copy()
		pConsv2.depends   = pBedRandom, pConsvPerm
		pConsv2.args.pval = True
		PyPPL(config).start(pBedRandom, pConsvPerm).run()
		procOK(pConsv2, 'consv.bed', self)


		
if __name__ == '__main__':
	unittest.main(failfast=True)