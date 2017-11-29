import unittest
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioaggrs.tfbs import aTfbsTfP, aTfbsP, aTfbsTfR, aTfbsR, aTfbsPC, aTfbsRC, aTfbsTfPC, aTfbsTfRC

class testTFBS (unittest.TestCase):

	def testTFBSOnPromotersByTF(self):
		aTfbsTfP.input = [[getfile('tfs.txt')], [getfile('genes.txt')]]
		aTfbsTfP.args.pval = 1e-3

		PyPPL().start(aTfbsTfP).run()
		procOK(aTfbsTfP.pMotifScan, 'aTfbsTfP.bed', self)

	def testTFBSOnPromoters(self):
		aTfbsP.input = [[getfile('tf-motifs.txt')], [getfile('genes.txt')]]
		aTfbsP.args.pval = 1e-3

		PyPPL().start(aTfbsP).run()
		procOK(aTfbsP.pMotifScan, 'aTfbsP.bed', self)

	def testTFBSOnRegionsByTF(self):
		aTfbsTfR.input = [[getfile('tfs.txt')], [getfile('regions.bed')]]
		aTfbsTfR.args.pval = 1e-3

		PyPPL().start(aTfbsTfR).run()
		procOK(aTfbsTfR.pMotifScan, 'aTfbsTfR.bed', self)

	def testTFBSOnRegions(self):
		aTfbsR.input = [[getfile('tf-motifs.txt')], [getfile('regions.bed')]]
		aTfbsR.args.pval = 1e-3

		PyPPL().start(aTfbsR).run()
		procOK(aTfbsR.pMotifScan, 'aTfbsR.bed', self)

	def testTFBSOnPromotersConsv(self):
		aTfbsPC.input = [[getfile('tf-motifs.txt')], [getfile('genes.txt')]]
		aTfbsPC.args.pval = 1e-3

		PyPPL().start(aTfbsPC).run()
		procOK(aTfbsPC.pConsv, 'aTfbsPC.bed', self)

	def testTFBSOnRegionsConsv(self):
		aTfbsRC.input = [[getfile('tf-motifs.txt')], [getfile('regions.bed')]]
		aTfbsRC.args.pval = 1e-3
		aTfbsRC.args.cpval = 0.1

		PyPPL().start(aTfbsRC).run()
		procOK(aTfbsRC.pConsv, 'aTfbsRC.bed', self)

	def testTFBSOnPromotersConsvByTF(self):
		aTfbsTfPC.input = [[getfile('tfs.txt')], [getfile('genes.txt')]]
		aTfbsTfPC.args.pval = 1e-3

		PyPPL().start(aTfbsTfPC).run()
		procOK(aTfbsTfPC.pConsv, 'aTfbsTfPC.bed', self)

	def testTFBSOnRegionsConsvByTF(self):
		aTfbsTfRC.input = [[getfile('tfs.txt')], [getfile('regions.bed')]]
		aTfbsTfRC.args.pval = 1e-3

		PyPPL().start(aTfbsTfRC).run()
		procOK(aTfbsTfRC.pConsv, 'aTfbsTfRC.bed', self)

if __name__ == '__main__':
	unittest.main(failfast = True)