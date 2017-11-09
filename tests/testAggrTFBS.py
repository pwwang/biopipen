import addPath, unittest

from os import path
from pyppl import PyPPL, params
from bioprocs.common import pStr2File
from bioaggrs.tfbs import aTFBSOnPromoters, aTFBSOnRegions, aTFBSOnPromotersConsv, aTFBSOnRegionsConsv, aTFBSOnPromotersByTF, aTFBSOnRegionsByTF
from bioprocs.bed import pBedRandom

mfile = params.tfmotifs.value
tffile = path.join(params.tmpdir.value, 'tffile.txt')
with open(tffile, 'w') as fout:
	fout.write("TEAD1_full_2	TEAD1\n")
	fout.write("M1916_1.02	TEAD1\n")
	fout.write("M4047_1.02	TEAD1\n")
	fout.write("M4049_1.02	TEAD1\n")
	fout.write("M4050_1.02	TEAD1\n")
	fout.write("M5904_1.02	TEAD1\n")
	fout.write("M5905_1.02	TEAD1\n")
	fout.write("TEAD3_DBD_1	TEAD3\n")
	fout.write("TEAD3_DBD_2	TEAD3\n")
	fout.write("M4052_1.02	TEAD3\n")
	fout.write("M5906_1.02	TEAD3\n")
	fout.write("M5907_1.02	TEAD3\n")
	fout.write("M6508_1.02	TEAD3\n")

class testAggrTFBS (unittest.TestCase):

	@unittest.skipIf(not path.isfile(mfile), 'Motif file not exists.')
	def testTFBSOnPromotersByTF(self):

		pStr2File5 = pStr2File.copy()
		pStr2File5.input = ['TEAD1, TEAD3']
		pStr2File6 = pStr2File.copy()
		pStr2File6.input = ['CXCL1,EMMPRIN']
		aTFBSOnPromotersByTF.depends2 = pStr2File5, pStr2File6
		aTFBSOnPromotersByTF.args.pval = 1e-3

		PyPPL().start(pStr2File5, pStr2File6, aTFBSOnPromotersByTF).run()

	@unittest.skipIf(not path.isfile(mfile), 'Motif file not exists.')
	def testTFBSOnPromoters(self):
		pStr2File.input = ['CXCL1,EMMPRIN']

		aTFBSOnPromoters.pPromoters.depends = pStr2File
		aTFBSOnPromoters.input = [[tffile]]
		aTFBSOnPromoters.args.pval = 1e-3

		PyPPL().start(pStr2File, aTFBSOnPromoters).run()

	@unittest.skipIf(not path.isfile(mfile), 'Motif file not exists.')
	def testTFBSOnRegionsByTF(self):
		pStr2File7 = pStr2File.copy()
		pStr2File7.input = ['TEAD1, TEAD3']
		pStr2File8 = pStr2File.copy()
		pStr2File8.input = ['chr19\t568135\t570239, chr9\t570035\t572139']

		aTFBSOnRegionsByTF.depends2 = pStr2File7, pStr2File8
		aTFBSOnRegionsByTF.args.pval = 1e-3

		PyPPL().start(pStr2File7, pStr2File8, aTFBSOnRegionsByTF).run()

	@unittest.skipIf(not path.isfile(mfile), 'Motif file not exists.')
	def testTFBSOnRegions(self):
		pStr2File2 = pStr2File.copy()

		pStr2File2.input = ['chr19\t568135\t570239, chr9\t570035\t572139']

		aTFBSOnRegions.pBedGetfasta.depends = pStr2File2
		aTFBSOnRegions.input = [[tffile]]
		aTFBSOnRegions.args.pval = 1e-3

		PyPPL().start(pStr2File2, aTFBSOnRegions).run()

	@unittest.skipIf(not path.isfile(mfile), 'Motif file not exists.')
	def testTFBSOnPromotersConsv(self):
		pStr2File3 = pStr2File.copy()
		pStr2File3.input = ['CXCL1,EMMPRIN']

		aTFBSOnPromotersConsv.pPromoters.depends = pStr2File3
		aTFBSOnPromotersConsv.input = [[tffile]]
		aTFBSOnPromotersConsv.args.pval = 1e-3

		PyPPL().start(pStr2File3, aTFBSOnPromotersConsv).run()

	@unittest.skipIf(not path.isfile(mfile), 'Motif file not exists.')
	def testTFBSOnPromotersConsv(self):
		pStr2File4 = pStr2File.copy()
		pStr2File4.input = ['chr19\t568135\t570239, chr9\t570035\t572139']

		aTFBSOnRegionsConsv.pBedGetfasta.depends = pStr2File4
		aTFBSOnRegionsConsv.input = [[tffile]]
		aTFBSOnRegionsConsv.args.pval = 1e-3
		aTFBSOnRegionsConsv.args.threspval = 0.1

		PyPPL().start(pStr2File4, aTFBSOnRegionsConsv).run()

if __name__ == '__main__':
	unittest.main()