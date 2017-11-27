import unittest
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.gene import pGeneNameNorm, pGeneTss, pGenePromoters

class TestGene (unittest.TestCase):
	def testGeneNameNorm(self):
		pGeneNameNorm.input     = [getfile('genes.txt')]
		pGeneNameNorm.errhow    = 'terminate'
		pGeneNameNorm.args.skip = 1
		pGeneNameNorm.args.col  = 1
		PyPPL(config).start(pGeneNameNorm).run()
		procOK(pGeneNameNorm, 'gene-namenorm.txt', self)

	def testGeneTss(self):
		pGeneTss.input = [getfile('genes.txt')]
		pGeneTss.errhow    = 'terminate'
		pGeneTss.args.skip = 1
		pGeneTss.args.col  = 1
		PyPPL(config).start(pGeneTss).run()
		procOK(pGeneTss, 'gene-tss.txt', self)

	def testpGenePromoters(self):
		pGenePromoters.input = [getfile('genes.txt')]
		pGenePromoters.errhow    = 'terminate'
		pGenePromoters.args.skip = 1
		pGenePromoters.args.col  = 1
		PyPPL(config).start(pGenePromoters).run()
		procOK(pGenePromoters, 'gene-proms.txt', self)


if __name__ == '__main__':
	unittest.main(failfast=True)