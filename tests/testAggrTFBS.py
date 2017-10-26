import addPath, unittest

from os import path
from pyppl import PyPPL
from bioprocs.common import pStr2File
from bioaggrs.tfbs import aTFBSOnPromoters, aTFBSOnRegions, aTFBSOnPromotersConsv, aTFBSOnRegionsConsv
from bioprocs.bed import pBedRandom

mfile = '/home/m161047/junwenwang/shared/reference/motif_databases/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme'

class testAggrTFBS (unittest.TestCase):

	@unittest.skipIf(not path.isfile(mfile), 'Motif file not exists.')
	def testTFBSOnPromoters(self):
		pStr2File.input = ['CXCL1,EMMPRIN']

		aTFBSOnPromoters.pPromoters.depends = pStr2File
		aTFBSOnPromoters.input = [[mfile]]
		aTFBSOnPromoters.args.mids = 'TEAD1_HUMAN.H10MO.D,TEAD3_HUMAN.H10MO.D'
		aTFBSOnPromoters.args.pval = 1e-3

		PyPPL().start(pStr2File, aTFBSOnPromoters).run()

	@unittest.skipIf(not path.isfile(mfile), 'Motif file not exists.')
	def testTFBSOnRegions(self):
		pStr2File2 = pStr2File.copy()

		pStr2File2.input = ['chr19\t568135\t570239, chr9\t570035\t572139']

		aTFBSOnRegions.pBedGetfasta.depends = pStr2File2
		aTFBSOnRegions.input = [[mfile]]
		aTFBSOnRegions.args.mids = 'TEAD1_HUMAN.H10MO.D,TEAD3_HUMAN.H10MO.D'
		aTFBSOnRegions.args.pval = 1e-3

		PyPPL().start(pStr2File2, aTFBSOnRegions).run()

	@unittest.skipIf(not path.isfile(mfile), 'Motif file not exists.')
	def testTFBSOnPromotersConsv(self):
		pStr2File3 = pStr2File.copy()
		pStr2File3.input = ['CXCL1,EMMPRIN']

		aTFBSOnPromotersConsv.pPromoters.depends = pStr2File3
		aTFBSOnPromotersConsv.input = [[mfile]]
		aTFBSOnPromotersConsv.args.mids = 'TEAD1_HUMAN.H10MO.D,TEAD3_HUMAN.H10MO.D'
		aTFBSOnPromotersConsv.args.pval = 1e-3

		PyPPL().start(pStr2File3, aTFBSOnPromotersConsv).run()

	@unittest.skipIf(not path.isfile(mfile), 'Motif file not exists.')
	def testTFBSOnPromotersConsv(self):
		pStr2File4 = pStr2File.copy()
		pStr2File4.input = ['chr19\t568135\t570239, chr9\t570035\t572139']

		aTFBSOnRegionsConsv.pBedGetfasta.depends = pStr2File4
		aTFBSOnRegionsConsv.input = [[mfile]]
		aTFBSOnRegionsConsv.args.mids = 'TEAD1_HUMAN.H10MO.D,TEAD3_HUMAN.H10MO.D'
		aTFBSOnRegionsConsv.args.pval = 1e-3
		aTFBSOnRegionsConsv.args.threspval = 0.1

		PyPPL().start(pStr2File4, aTFBSOnRegionsConsv).run()

if __name__ == '__main__':
	unittest.main()