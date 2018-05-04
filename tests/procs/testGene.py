import helpers, testly
from os import path
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.gene import pGeneNameNorm, pGeneTss, pGenePromoters, pGeneBody

class TestGene (helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestGene')
	def testGeneNameNorm(self):
		pGeneNameNorm.input             = [path.join(self.indir, 'genes.txt')]
		pGeneNameNorm.errhow            = 'terminate'
		pGeneNameNorm.args.inopts.skip  = 1
		pGeneNameNorm.args.outopts.head = False
		pGeneNameNorm.args.genecol      = 'COL2'
		PyPPL(config).start(pGeneNameNorm).run()
		self.assertFileEqual(pGeneNameNorm.channel.get(), path.join(self.outdir, 'gene-namenorm.txt'))

	def testGeneTss(self):
		pGeneTss.input              = [getfile('genes.txt')]
		pGeneTss.errhow             = 'terminate'
		pGeneTss.echo   = 0
		pGeneTss.args.inopts.skip   = 1
		pGeneTss.args.genecol       = 'COL2'
		pGeneTss.args.outopts.query = True
		PyPPL(config).start(pGeneTss).run()
		self.assertFileCountEqual(pGeneTss.channel.get(), path.join(self.outdir, 'gene-tss.txt'))

	def testpGenePromoters(self):
		pGenePromoters.input              = [getfile('genes.txt')]
		pGenePromoters.errhow             = 'terminate'
		pGenePromoters.args.inopts.skip   = 1
		pGenePromoters.args.genecol       = 'COL2'
		pGenePromoters.args.outopts.query = True
		PyPPL(config).start(pGenePromoters).run()
		self.assertFileCountEqual(pGenePromoters.channel.get(), path.join(self.outdir, 'gene-proms.txt'))

	def testGeneTss(self):
		pGeneBody.input              = [getfile('genes.txt')]
		pGeneBody.errhow             = 'terminate'
		pGeneBody.echo   = 0
		pGeneBody.args.inopts.skip   = 1
		pGeneBody.args.genecol       = 'COL2'
		pGeneBody.args.outopts.query = True
		PyPPL(config).start(pGeneBody).run()
		self.assertFileCountEqual(pGeneBody.channel.get(), path.join(self.outdir, 'gene-body.bedx'))


if __name__ == '__main__':
	testly.main(failfast=True)
