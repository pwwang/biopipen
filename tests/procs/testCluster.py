import unittest
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.cluster import pDist2Coords, pCluster, pMCluster, pAPCluster, pHCluster

class testCluster (unittest.TestCase):
	
	def testDist2Coords (self):
		pDist2Coords1       = pDist2Coords.copy()
		pDist2Coords1.input = [getfile('dist2coords.in.full.txt')]
		PyPPL(config).start(pDist2Coords1).run()
		procOK(pDist2Coords1, 'dist2coords.out.txt', self)

		pDist2Coords2            = pDist2Coords.copy()
		pDist2Coords2.input      = [getfile('dist2coords.in.upper.txt')]
		pDist2Coords2.args.infmt = 'upper'
		PyPPL(config).start(pDist2Coords2).run()
		procOK(pDist2Coords2, 'dist2coords.out.txt', self)

		pDist2Coords3            = pDist2Coords.copy()
		pDist2Coords3.input      = [getfile('dist2coords.in.lower.txt')]
		pDist2Coords3.args.infmt = 'lower'
		PyPPL(config).start(pDist2Coords3).run()
		procOK(pDist2Coords3, 'dist2coords.out.txt', self)

		pDist2Coords4            = pDist2Coords.copy()
		pDist2Coords4.input      = [getfile('dist2coords.in.pair.txt')]
		pDist2Coords4.args.infmt = 'pair'
		PyPPL(config).start(pDist2Coords4).run()
		procOK(pDist2Coords4, 'dist2coords.out.txt', self, order=False)

	def testCluster(self):
		pCluster.input     = [getfile('cluster.txt')]
		pCluster.args.minc = 4
		pCluster.args.maxc = 6
		PyPPL(config).start(pCluster).run()
		procOK(pCluster, 'cluster.txt', self)
	
	def testMCluster(self):
		pMCluster.input     = [getfile('cluster.txt')]
		pMCluster.args.minc = 4
		pMCluster.args.maxc = 6
		PyPPL(config).start(pMCluster).run()
		procOK(pMCluster, 'cluster.txt', self)
	
	def testAPCluster(self):
		pAPCluster.input     = [getfile('cluster.txt')]
		PyPPL(config).start(pAPCluster).run()
		procOK(pAPCluster, 'apcluster.txt', self, order=False)

	def testHCluster(self):
		pHCluster.input     = [getfile('cluster.txt')]
		PyPPL(config).start(pHCluster).run()
		procOK(pHCluster, 'hclust-merge.txt', self)

if __name__ == '__main__':
	unittest.main()
		