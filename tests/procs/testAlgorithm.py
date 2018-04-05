import helpers, unittest
from os import path
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.algorithm import pRWR, pAR

class testAlgorithm (helpers.TestCase):
	
	def dataProvider_testRWR(self, testdir, indir, outdir):
		wfile = path.join(indir, 'test.rwr.w.txt')
		efile = path.join(indir, 'test.rwr.e.txt')
		outfile = path.join(outdir, 'test.rwr.txt')
		yield 't1', wfile, efile, outfile
	
	def testRWR (self, tag, wfile, efile, outfile):
		self.maxDiff = None
		pRWRTest = pRWR.copy(tag = tag)
		pRWRTest.input = wfile, efile
		PyPPL(config).start(pRWRTest).run()
		self.assertFileEqual(pRWRTest.channel.get(), outfile)
		
	def dataProvider_testAR(self, testdir, indir, outdir):
		dfile  = path.join(indir, 'ar.d.txt.gz')
		ptfile = path.join(indir, 'ar.pt.txt.gz')
		yfile  = path.join(indir, 'ar.y.txt.gz')
		wfile  = path.join(outdir, 'ar.w.txt')
		args   = {'parallel': True, 'svdP': 25}
		yield 't1', dfile, ptfile, yfile, args, wfile
		args1  = {'parallel': True, 'svdP': 25, 'method': 'admm'}
		#ADMM take too long
		#yield 't2', dfile, ptfile, yfile, args1, wfile

	def testAR(self, tag, dfile, ptfile, yfile, args, wfile):
		pARTest = pAR.copy(tag = tag)
		pARTest.input = dfile, ptfile, yfile
		pARTest.args.update(args)
		PyPPL(config).start(pARTest).run()
		self.assertFileEqual(pARTest.channel.get(), wfile)
	
if __name__ == '__main__':
	unittest.main()
		