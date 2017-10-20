import addPath, unittest

from pyppl import PyPPL
from bioprocs.common import pStr2File
from bioaggrs.tfbs import aTFBSOnPromoters, aTFBSOnRegions

class testAggrTFBS (unittest.TestCase):

	def testTFBSOnPromoters(self):
		pStr2File.input = ['CXCL1,EMMPRIN']

		aTFBSOnPromoters.pPromoters.depends = pStr2File
		aTFBSOnPromoters.input = [['/home/m161047/junwenwang/shared/reference/motif_databases/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme']]
		aTFBSOnPromoters.args.mids = 'TEAD1_HUMAN.H10MO.D,TEAD3_HUMAN.H10MO.D'
		aTFBSOnPromoters.args.pval = 1e-3

		PyPPL().start(pStr2File, aTFBSOnPromoters).run()

	def testTFBSOnRegions(self):
		pStr2File2 = pStr2File.copy()

		pStr2File2.input = ['chr19\t568135\t570239, chr9\t570035\t572139']

		aTFBSOnRegions.pBedGetfasta.depends = pStr2File2
		aTFBSOnRegions.input = [['/home/m161047/junwenwang/shared/reference/motif_databases/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme']]
		aTFBSOnRegions.args.mids = 'TEAD1_HUMAN.H10MO.D,TEAD3_HUMAN.H10MO.D'
		aTFBSOnRegions.args.pval = 1e-3

		PyPPL().start(pStr2File2, aTFBSOnRegions).run()

if __name__ == '__main__':
	unittest.main()