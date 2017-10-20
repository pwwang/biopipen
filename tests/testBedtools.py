import addPath, unittest
from pyppl import PyPPL
from os import path
from tempfile import gettempdir
from bioprocs.bedtools import pBedGetfasta
from bioprocs.common import pStr2File

def readFile(file):
	with open(file) as f:
		return [line.strip() for line in f if line.strip()]

def readProc(proc):
	return readFile(proc.channel.get())

tmpdir = gettempdir()
class testBedtools (unittest.TestCase):
	
	def testGetfasta(self):
		pStr2File1 = pStr2File.copy()
		pStr2File1.input = ['chr1	5	10']
		pStr2File2 = pStr2File.copy()
		pStr2File2.input = ['>chr1,AAAAAAAACCCCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG']

		pBedGetfasta.depends = pStr2File1, pStr2File2
		pBedGetfasta.callfront = lambda p: p.args.update({'ref': pStr2File2.channel.get()})

		PyPPL().start(pStr2File1, pStr2File2).run()
		self.assertEqual(readProc(pBedGetfasta), ['>::chr1:5-10', 'AAACC'])
	'''
	def testClosest (self):
		afile  = os.path.join("testfiles", "bedtools", "test.closest.a.bed")
		bfile1 = os.path.join("testfiles", "bedtools", "test.closest.b1.bed")
		bfile2 = os.path.join("testfiles", "bedtools", "test.closest.b2.bed")
		outfile = os.path.join("testfiles", "bedtools", "test.closest.ret")
		pClosest.input = {pClosest.input: [(afile, [bfile1, bfile2])]}
		pClosest.run ()
		self.assertEqual (readFile(outfile), readProc(pClosest))
		
	def testFlank (self):
		infile  = os.path.join("testfiles", "bedtools", "test.flank.bed")
		gfile   = os.path.join("testfiles", "bedtools", "test.flank.genome")
		outfile = os.path.join("testfiles", "bedtools", "test.flank.ret")
		pFlank.input = {pFlank.input: [(infile, gfile)]}
		pFlank.args.update ({"params": "-b 5"})
		pFlank.run ()
		self.assertEqual (readFile(outfile), readProc(pFlank))
		
	def testIntersect (self):
		afile  = os.path.join("testfiles", "bedtools", "test.intersect.a.bed")
		bfile  = os.path.join("testfiles", "bedtools", "test.intersect.b.bed")
		outfile = os.path.join("testfiles", "bedtools", "test.intersect.ret")
		pIntersect.input = {pIntersect.input: [(afile, bfile)]}
		pIntersect.run ()
		self.assertEqual (readFile(outfile), readProc(pIntersect))
		
	def testMakewindows (self):
		infile  = os.path.join("testfiles", "bedtools", "test.makewindows.bed")
		outfile = os.path.join("testfiles", "bedtools", "test.makewindows.ret")
		pMakewindows.input = {pMakewindows.input: [infile]}
		pMakewindows.args.update ({"params": "-n 5", "informat": "bed"})
		pMakewindows.run ()
		self.assertEqual (readFile(outfile), readProc(pMakewindows))
		
	def testMerge (self):
		infile  = os.path.join("testfiles", "bedtools", "test.merge.bed")
		outfile = os.path.join("testfiles", "bedtools", "test.merge.ret")
		pMerge.input = {pMerge.input: [infile]}
		pMerge.run ()
		self.assertEqual (readFile(outfile), readProc(pMerge))
		
	def testMultiinter (self):
		infilea  = os.path.join("testfiles", "bedtools", "test.multiinter.a.bed")
		infileb  = os.path.join("testfiles", "bedtools", "test.multiinter.b.bed")
		infilec  = os.path.join("testfiles", "bedtools", "test.multiinter.c.bed")
		outfile = os.path.join("testfiles", "bedtools", "test.multiinter.ret")
		pMultiinter.input = {pMultiinter.input: [([infilea, infileb, infilec], )]}
		pMultiinter.run ()
		self.assertEqual (readFile(outfile), readProc(pMultiinter))
		
	def testRandom (self):
		gfile  = os.path.join("testfiles", "bedtools", "test.random.genome")
		outfile = os.path.join("testfiles", "bedtools", "test.random.ret")
		pRandom.input = {pRandom.input: [gfile]}
		pRandom.args.update({ "params": "-n 10 -seed 0" })
		pRandom.run ()
		self.assertEqual (readFile(outfile), readProc(pRandom))
		
	def testShift (self):
		infile = os.path.join("testfiles", "bedtools", "test.shift.bed")
		gfile  = os.path.join("testfiles", "bedtools", "test.shift.genome")
		outfile = os.path.join("testfiles", "bedtools", "test.shift.ret")
		pShift.input = {pShift.input: [(infile, gfile)]}
		pShift.args.update({ "params": "-s 5" })
		pShift.run ()
		self.assertEqual (readFile(outfile), readProc(pShift))
		
	def testShuffle (self):
		infile = os.path.join("testfiles", "bedtools", "test.shuffle.bed")
		gfile  = os.path.join("testfiles", "bedtools", "test.shuffle.genome")
		outfile = os.path.join("testfiles", "bedtools", "test.shuffle.ret")
		pShuffle.input = {pShuffle.input: [(infile, gfile)]}
		pShuffle.args.update({ "params": "-seed 0" })
		pShuffle.run ()
		self.assertEqual (readFile(outfile), readProc(pShuffle))
		
	def testSubtract (self):
		afile = os.path.join("testfiles", "bedtools", "test.subtract.a.bed")
		bfile = os.path.join("testfiles", "bedtools", "test.subtract.b.bed")
		outfile = os.path.join("testfiles", "bedtools", "test.subtract.ret")
		pSubtract.input = {pSubtract.input: [(afile, bfile)]}
		pSubtract.args.update({ "params": "-seed 0" })
		pSubtract.run ()
		self.assertEqual (readFile(outfile), readProc(pSubtract))
		
	def testWindow (self):
		afile = os.path.join("testfiles", "bedtools", "test.window.a.bed")
		bfile = os.path.join("testfiles", "bedtools", "test.window.b.bed")
		outfile = os.path.join("testfiles", "bedtools", "test.window.ret")
		pWindow.input = {pWindow.input: [(afile, bfile)]}
		pWindow.run ()
		self.assertEqual (readFile(outfile), readProc(pWindow))
	'''
if __name__ == '__main__':
	unittest.main()
		