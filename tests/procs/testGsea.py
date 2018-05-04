import unittest, helpers
from os import path
from pyppl import PyPPL
from helpers import getfile, procOK, config, testdirs
from bioprocs.gsea import pTargetEnrichr, pEnrichr, pExpmat2Gct, pSSGSEA, pGSEA, pSampleinfo2Cls, pGMT2Mat


class TestGsea (helpers.TestCase):

	testdir, indir, outdir = testdirs('TestGsea')

	def testpTargetEnrichrCol2 (self):
		pTargetEnrichrCol2        = pTargetEnrichr.copy()
		pTargetEnrichrCol2.input  = [getfile('tenrich-col2.txt')]
		pTargetEnrichrCol2.errhow = 'terminate'
		PyPPL(config).start(pTargetEnrichrCol2).run()
		self.assertFileEqual(path.join(pTargetEnrichrCol2.channel.get(), 'tenrich-col2-KEGG_2016.txt'), path.join(self.outdir, 'tenrich-col2-KEGG_2016.txt'))

	def testpTargetEnrichrCol3 (self):
		pTargetEnrichrCol3        = pTargetEnrichr.copy()
		pTargetEnrichrCol3.input  = [getfile('tenrich-col3.txt')]
		pTargetEnrichrCol3.errhow = 'terminate'
		PyPPL(config).start(pTargetEnrichrCol3).run()
		self.assertFileEqual(path.join(pTargetEnrichrCol3.channel.get(), 'tenrich-col3-KEGG_2016.txt'), path.join(self.outdir, 'tenrich-col3-KEGG_2016.txt'))

	def testpTargetEnrichrCol4 (self):
		pTargetEnrichrCol4        = pTargetEnrichr.copy()
		pTargetEnrichrCol4.input  = [getfile('tenrich-col4.txt')]
		pTargetEnrichrCol4.errhow = 'terminate'
		PyPPL(config).start(pTargetEnrichrCol4).run()
		self.assertFileEqual(path.join(pTargetEnrichrCol4.channel.get(), 'tenrich-col4-KEGG_2016.txt'), path.join(self.outdir, 'tenrich-col4-KEGG_2016.txt'))

	def testpTargetEnrichrCol5 (self):
		pTargetEnrichrCol5        = pTargetEnrichr.copy()
		pTargetEnrichrCol5.input  = [getfile('tenrich-col5.txt')]
		pTargetEnrichrCol5.errhow = 'terminate'
		PyPPL(config).start(pTargetEnrichrCol5).run()
		self.assertFileEqual(path.join(pTargetEnrichrCol5.channel.get(), 'tenrich-col5-KEGG_2016.txt'), path.join(self.outdir, 'tenrich-col5-KEGG_2016.txt'))

	def testpEnrichr(self):
		self.maxDiff = None
		pEnrichr.input     = [getfile('enrichr.txt')]
		pEnrichr.args.norm = True
		pEnrichr.errhow    = 'terminate'
		PyPPL(config).start(pEnrichr).run()
		self.assertFileEqual(path.join(pEnrichr.channel.get(), 'enrichr-KEGG_2016.txt'), path.join(self.outdir, 'enrichr-KEGG_2016.txt'))

	def testpExpmat2Gct(self):
		pExpmat2GctCopy   = pExpmat2Gct.copy()
		pExpmat2GctCopy.input = [getfile('expmat.txt')]
		PyPPL(config).start(pExpmat2GctCopy).run()
		procOK(pExpmat2GctCopy, 'expmat.txt', self)

	def testpSSGSEA(self):
		self.maxDiff = None
		pExpmat2GctCopy2       = pExpmat2Gct.copy()
		pExpmat2GctCopy2.input = [getfile('ssgct.txt')]
		pSSGSEA.depends        = pExpmat2GctCopy2
		pSSGSEA.input          = lambda ch: ch.cbind(getfile('ssgmt.txt'))
		pSSGSEA.args.seed      = 8525
		PyPPL(config).start(pExpmat2GctCopy2).run()
		self.assertFileEqual(path.join(pSSGSEA.channel.get(), 'report.txt'), path.join(self.outdir, 'ssgsea-report.txt'))

	def testpGSEA(self):
		pSampleinfo2Cls.input  = [getfile('saminfo.txt')]
		pExpmat2GctCopy3       = pExpmat2Gct.copy()
		pExpmat2GctCopy3.input = [getfile('gct.txt')]
		pGSEA.depends          = pExpmat2GctCopy3, pSampleinfo2Cls
		pGSEA.args.nperm       = 1000
		pGSEA.args.nthread     = 4
		pGSEA.input            = lambda ch,   ch2: ch.cbind(ch2, getfile('gmt.txt'))
		PyPPL(config).start(pExpmat2GctCopy3, pSampleinfo2Cls).run()
		procOK(pGSEA, 'gsea.txt', self)

	def testpGMT2Mat(self):
		pGMT2Mat.input = [getfile('gmt2mat.txt')]
		PyPPL(config).start(pGMT2Mat).run()
		procOK(pGMT2Mat, 'gmt2mat.txt', self)

if __name__ == '__main__':
	unittest.main(failfast=True)
