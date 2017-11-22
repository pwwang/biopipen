import unittest
from os import path
from pyppl import PyPPL
from bioprocs.tfbs import pMotifScan
from helpers import getfile, procOK, config

class testTfbs (unittest.TestCase):
	
	def testMotifScanOnlyMotif(self):
		name1             = 'motifs.txt'
		name2             = 'motifs.fa'
		pMotifScan1       = pMotifScan.copy()
		pMotifScan1.input = (getfile(name1), getfile(name2))
		PyPPL(config).start(pMotifScan1).run()
		procOK(pMotifScan1, 'motifs-motifs.bed', self)
	
	def testMotifScanTF(self):
		name1             = 'tf-motifs.txt'
		name2             = 'motifs.fa'
		pMotifScan2       = pMotifScan.copy()
		pMotifScan2.input = (getfile(name1), getfile(name2))
		PyPPL(config).start(pMotifScan2).run()
		procOK(pMotifScan2, 'tf-motifs-motifs.bed', self)
	
if __name__ == '__main__':
	unittest.main()
		