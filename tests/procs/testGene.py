import helpers, testly
from os import path
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.gene import pGeneNameNorm, pGeneTss, pGenePromoters

class TestGene (helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestGene')
	def testGeneNameNorm(self):
		pGeneNameNorm.input     = [path.join(self.indir, 'genes.txt')]
		pGeneNameNorm.errhow    = 'terminate'
		pGeneNameNorm.args.inopts.skip = 1
		pGeneNameNorm.args.outopts.head = False
		pGeneNameNorm.args.genecol = 'COL2'
		PyPPL(config).start(pGeneNameNorm).run()
		self.assertFileEqual(pGeneNameNorm.channel.get(), path.join(self.outdir, 'gene-namenorm.txt'))

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
	testly.main(failfast=True)
