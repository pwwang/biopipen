import sys, os, json

sys.path.insert(0, "/home/m161047/tools/pyppl")
sys.path.insert(0, "/home/m161047/tools/bioprocs")

import pyppl, unittest
from bioprocs.algorithm import *

def readFile(file):
	return [line.strip() for line in open(file)]

def readProc(proc):
	return [line.strip() for line in open(proc.output['outfile'][0])]

class testAlgorithm (unittest.TestCase):
	
	def testRWR (self):
		self.maxDiff = None
		wfile   = os.path.join ("testfiles", "algorithm", "test.rwr.w.txt")
		efile   = os.path.join ("testfiles", "algorithm", "test.rwr.e.txt")
		outfile = os.path.join ("testfiles", "algorithm", "test.rwr.ret")
		pRWR.input = {pRWR.input: [(wfile, efile)]}
		pRWR.args.update({
			'Wformat': 'txt',
			'Eformat': 'txt',
			'Rformat': 'txt',
			'normW': True,
			'normE': True
		})
		pRWR.run()
		self.assertEqual (readFile(outfile), readProc(pRWR))
		
if __name__ == '__main__':
	unittest.main()
		