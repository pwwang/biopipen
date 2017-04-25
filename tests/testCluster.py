import sys, os, json

sys.path.insert(0, "/home/m161047/tools/pyppl")
sys.path.insert(0, "/home/m161047/tools/bioprocs")

import pyppl, unittest
from bioprocs.cluster import *

def readFile(file):
	clusters = {}
	with open (file) as f:
		for line in f:
			(i, c) = line.split("\t")
			if clusters.has_key (c):
				clusters[c].append (i)
			else:
				clusters[c] = [i]
	return sorted(clusters.values())

class testCluster (unittest.TestCase):
	
	def testDist2Coords (self):
		infile1 = "testfiles/cluster/dist2coords.in.full.txt"
		infile2 = "testfiles/cluster/dist2coords.in.triangle.txt"
		infile3 = "testfiles/cluster/dist2coords.in.pair.txt"
		outfile = "testfiles/cluster/dist2coords.out.txt"
		
		inkeys  = pDist2Coords.input
		pDist2Coords.input  = {inkeys: [infile1]}
		pDist2Coords.args['informat'] = 'full'
		pDist2Coords.run()
		self.assertEqual (readFile(outfile), readFile(pDist2Coords.output['outfile'][0]))
		
		pDist2Coords.input  = {inkeys: [infile2]}
		pDist2Coords.args['informat'] = 'triangle'
		pDist2Coords.run()
		self.assertEqual (readFile(outfile), readFile(pDist2Coords.output['outfile'][0]))
		
		pDist2Coords.input  = {inkeys: [infile3]}
		pDist2Coords.args['informat'] = 'pair'
		pDist2Coords.run()
		self.assertEqual (readFile(outfile), readFile(pDist2Coords.output['outfile'][0]))
		
	def testDecideK (self):
		infile   = "testfiles/cluster/kmeans.input.txt"
		pDecideK.input = {pDecideK.input: [infile]}
		pDecideK.args['method'] = 'elbow'
		pDecideK.tag = pDecideK.args['method']
		pDecideK.run ()
		self.assertEqual (open(pDecideK.output['kfile'][0]).read().strip(), "7")
		
		pDecideK.args['method'] = "pamk"
		pDecideK.tag = pDecideK.args['method']
		pDecideK.run ()
		self.assertEqual (open(pDecideK.output['kfile'][0]).read().strip(), "12")
		
		pDecideK.args['method'] = "calinski"
		pDecideK.tag = pDecideK.args['method']
		pDecideK.run ()
		self.assertEqual (open(pDecideK.output['kfile'][0]).read().strip(), "13")
		
		pDecideK.args['method'] = "mclust"
		pDecideK.tag = pDecideK.args['method']
		pDecideK.run ()
		self.assertEqual (open(pDecideK.output['kfile'][0]).read().strip(), "7")
		
		pDecideK.args['method'] = "ap"
		pDecideK.tag = pDecideK.args['method']
		pDecideK.run ()
		self.assertEqual (open(pDecideK.output['kfile'][0]).read().strip(), "7")
		
		pDecideK.args['method'] = "gap"
		pDecideK.tag = pDecideK.args['method']
		pDecideK.run ()
		self.assertEqual (open(pDecideK.output['kfile'][0]).read().strip(), "7")
		
		pDecideK.args['method'] = "nbclust"
		pDecideK.tag = pDecideK.args['method']
		pDecideK.run ()
		self.assertEqual (open(pDecideK.output['kfile'][0]).read().strip(), "6")
	
	def testKMeans (self):
		infile  = "testfiles/cluster/kmeans.input.txt"
		outfile = "testfiles/cluster/kmeans.output.txt"
		
		pKMeans.input = {pKMeans.input: [(infile, 7)]}
		pKMeans.args['row.names'] = 1
		pKMeans.args['header']    = True
		pKMeans.run()
		self.assertEqual (readFile(outfile), readFile(pKMeans.output['outdir'][0] + '/cluster.txt'))
		
	def testPamk (self):
		infile  = "testfiles/cluster/kmeans.input.txt"
		outfile = "testfiles/cluster/pamk.out.txt"
		pPamk.input = {pPamk.input: [infile]}
		pPamk.args['row.names'] = 1
		pPamk.args['header']    = True
		pPamk.run()
		self.assertEqual (readFile(outfile), readFile(pPamk.output['outdir'][0] + '/cluster.txt'))
		
	def testClara (self):
		infile  = "testfiles/cluster/kmeans.input.txt"
		kfile   = "testfiles/cluster/clara.k.txt"
		outfile = "testfiles/cluster/clara.out.txt"
		pClara.input  = {pClara.input: [(infile, kfile)]}
		pClara.run()
		self.assertEqual (readFile(outfile), readFile(pClara.output['outdir'][0] + '/cluster.txt'))
		
	def testMClust (self):
		infile  = "testfiles/cluster/kmeans.input.txt"
		outfile = "testfiles/cluster/clara.out.txt"
		pMClust.input  = {pMClust.input: [(infile)]}
		pMClust.run()
		self.assertEqual (readFile(outfile), readFile(pMClust.output['outdir'][0] + '/cluster.txt'))
		
	def testAPCluster (self):
		infile  = "testfiles/cluster/kmeans.input.txt"
		outfile = "testfiles/cluster/clara.out.txt"
		pAPCluster.input  = {pAPCluster.input: [(infile)]}
		pAPCluster.run()
		self.assertEqual (readFile(outfile), readFile(pAPCluster.output['outdir'][0] + '/cluster.txt'))
		
	def testHClust (self):
		infile  = "testfiles/cluster/kmeans.input.txt"
		mergefile = "testfiles/cluster/hclust.merge.txt"
		orderfile = "testfiles/cluster/hclust.order.txt"
		pHClust.input  = {pHClust.input: [(infile)]}
		pHClust.args['gg'] = True
		pHClust.run()
		self.assertEqual (open(mergefile).read(), open(pHClust.output['outdir'][0] + '/hclust.merge.txt').read())
		self.assertEqual (open(orderfile).read(), open(pHClust.output['outdir'][0] + '/hclust.order.txt').read())
	
if __name__ == '__main__':
	unittest.main()
		