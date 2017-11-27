import unittest
from pyppl import PyPPL
from helpers import getfile, procOK, config
from pyppl import PyPPL
from bioprocs.snp import pSnp2Bed, pSnp2Avinput


class TestSnp(unittest.TestCase):
	def testpSnp2Bed(self):
		pSnp2Bed.input     = [getfile('snps.txt')]
		pSnp2Bed.errhow    = 'terminate'
		pSnp2Bed.args.skip = 0
		pSnp2Bed.args.col  = 0
		PyPPL(config).start(pSnp2Bed).run()
		procOK(pSnp2Bed, 'snps.bed', self)

	def testpSnp2Avinput(self):
		pSnp2Avinput.input = [getfile('snps.txt')]
		pSnp2Avinput.errhow    = 'terminate'
		pSnp2Avinput.args.skip = 0
		pSnp2Avinput.args.col  = 0
		PyPPL(config).start(pSnp2Avinput).run()
		procOK(pSnp2Avinput, 'snps.avinput', self)

if __name__ == '__main__':
	unittest.main(failfast=True)