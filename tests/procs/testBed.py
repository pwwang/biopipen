import unittest
from pyppl import PyPPL
from os import path
from bioprocs.bed import pBedGetfasta, pBedFlank, pBedRandom, pBedSort
from helpers import getfile, procOK, config

class testBedtools (unittest.TestCase):
	
	def testGetfasta(self):
		pBedGetfasta.input    = [getfile('bedgetfasta.bed')]
		pBedGetfasta.args.ref = getfile('bedgetfasta.fa')

		PyPPL(config).start(pBedGetfasta).run()
		procOK(pBedGetfasta, 'bedgetfasta.fa', self)
		
	def testFlank (self):
		pBedFlank.input         = [getfile('bedflank.bed')]
		pBedFlank.args.params.b = 5
		pBedFlank2              = pBedFlank.copy()
		pBedFlank2.input        = [getfile('bedflank.bed')]
		pBedFlank2.args.extend  = True
		PyPPL().start(pBedFlank, pBedFlank2).run()
		procOK(pBedFlank, 'bedflank1.bed', self)
		procOK(pBedFlank2, 'bedflank2.bed', self)

	def testRandom(self):
		pBedRandom.input     = (100, 10)
		pBedRandom.args.seed = 8525

		PyPPL(config).start(pBedRandom).run()
		procOK(pBedRandom, 'bedrandom.bed', self)

	def testGetfasta(self):
		pBedSort.input = [getfile('bedsort.bed')]

		pBedSortSort             = pBedSort.copy()
		pBedSortSort.args.tool   = 'sort'
		pBedSortSort.args.unique = False

		pBedSortSortU             = pBedSort.copy()
		pBedSortSortU.args.tool   = 'sort'
		pBedSortSortU.args.unique = True

		pBedSortBedOps             = pBedSort.copy()
		pBedSortBedOps.args.tool   = 'bedops'
		pBedSortBedOps.args.unique = False

		pBedSortBedOpsU             = pBedSort.copy()
		pBedSortBedOpsU.args.tool   = 'bedops'
		pBedSortBedOpsU.args.unique = True

		pBedSortBedtools             = pBedSort.copy()
		pBedSortBedtools.args.tool   = 'bedtools'
		pBedSortBedtools.args.unique = False

		pBedSortBedtoolsU             = pBedSort.copy()
		pBedSortBedtoolsU.args.tool   = 'bedtools'
		pBedSortBedtoolsU.args.unique = True

		PyPPL(config).start(pBedSortSort, pBedSortSortU, pBedSortBedOps, pBedSortBedOpsU, pBedSortBedtools, pBedSortBedtoolsU).run()
		procOK(pBedSortSort, 'bedsort.bed', self)
		procOK(pBedSortSortU, 'bedsort-unique.bed', self)
		procOK(pBedSortBedOps, 'bedsort.bed', self)
		procOK(pBedSortBedOpsU, 'bedsort-unique.bed', self)
		procOK(pBedSortBedtools, 'bedsort.bed', self)
		procOK(pBedSortBedtoolsU, 'bedsort-unique.bed', self)
	
	# TODO: 
	# - pBedClosest
	# - pBedIntersect
	# - pBedMakewindows
	# - pBedMerge
	# - pBedMultiinter
	# - pBedShift
	# - pBedShuffle
	# - pBedSubtract
	# - pBedWindow
	# - pBedGenomecov
	# - pBedCluster
if __name__ == '__main__':
	unittest.main()
		