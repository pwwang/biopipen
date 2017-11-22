import unittest
from os import path
from pyppl import PyPPL
from bioprocs.matrix import pMatrix, pCbind, pRbind, pCsplit, pRsplit, pTxtFilter, pTxtTransform, pSimRead
from helpers import getfile, procOK, config

class testMatrix (unittest.TestCase):
	
	def testMatrix(self):
		name               = 'matrix.txt'
		pMatrix1           = pMatrix.copy()
		pMatrix1.input     = [getfile(name)]
		pMatrix1.args.code = "mat = log(mat + 1, 2)"
		PyPPL(config).start(pMatrix1).run()
		procOK(pMatrix1, name, self)

	def testMatrixNoHead(self):
		name                 = 'matrix-nohead.txt'
		pMatrix2             = pMatrix.copy()
		pMatrix2.input       = [getfile(name)]
		pMatrix2.args.code   = "mat = log(mat + 1, 2)"
		pMatrix2.args.cnames = False
		PyPPL(config).start(pMatrix2).run()
		procOK(pMatrix2, name, self)

	def testMatrixNoRn(self):
		name                 = 'matrix-norn.txt'
		pMatrix3             = pMatrix.copy()
		pMatrix3.input       = [getfile(name)]
		pMatrix3.args.code   = "mat = log(mat + 1, 2)"
		pMatrix3.args.rnames = False
		PyPPL(config).start(pMatrix3).run()
		procOK(pMatrix3, name, self)

	def testCbind(self):
		name1 = 'matrix-cbind-1.txt'
		name2 = 'matrix-cbind-2.txt'
		name0 = 'matrix-cbind.txt'
		pCbind1 = pCbind.copy()
		pCbind1.input = [[getfile(name1), getfile(name2)]]
		PyPPL(config).start(pCbind1).run()
		procOK(pCbind1, name0, self)

	def testRbind(self):
		name1 = 'matrix-rbind-1.txt'
		name2 = 'matrix-rbind-2.txt'
		name0 = 'matrix-rbind.txt'
		pRbind1 = pRbind.copy()
		pRbind1.input = [[getfile(name1), getfile(name2)]]
		PyPPL(config).start(pRbind1).run()
		procOK(pRbind1, name0, self)

	def testCsplit(self):
		name1 = 'matrix-rbind-1.txt'
		pCsplit.input = [getfile(name1)]
		PyPPL(config).start(pCsplit).run()
		procOK(pCsplit, 'matrix-rbind-1.splits', self)

	def testRsplit(self):
		name1 = 'matrix-rbind-2.txt'
		pRsplit.input = [getfile(name1)]
		PyPPL(config).start(pRsplit).run()
		procOK(pRsplit, 'matrix-rbind-2.splits', self)

	def testTxtFilter(self):
		pTxtFilter.input        = [getfile('txtfilter.txt')]
		pTxtFilter.args.cols    = [0, 2, "V4"]
		pTxtFilter.args.rfilter = 'lambda row: float(row[1]) > 2 and float(row[2]) > 2'
		pTxtFilter.args.skip    = 2
		pTxtFilter.args.delimit = '|'
		PyPPL(config).start(pTxtFilter).run()
		procOK(pTxtFilter, 'txtfilter.txt', self)

	def testTxtTransform(self):
		pTxtTransform.input          = [getfile('txtfilter.txt')]
		pTxtTransform.args.cols      = [0, 2, "V4"]
		pTxtTransform.args.transform = 'lambda row: [str(float(r) + 1) if i == 1 else r for i, r in enumerate(row)]'
		pTxtTransform.args.skip      = 2
		pTxtTransform.args.delimit   = '|'
		PyPPL(config).start(pTxtTransform).run()
		procOK(pTxtTransform, 'txttransform.txt', self)

	def testSimRead(self):
		pSimRead.input        = [[getfile('simread1.txt.gz'), getfile('simread2.txt')]]
		pSimRead.args.skip    = [3]
		pSimRead.args.usehead = 0
		pSimRead.args.match   = 'lambda r1, r2: -1 if r1[0] == r2[1] else 0 if r1[0] < r2[1] else 1'
		pSimRead.args.do      = 'lambda r1, r2: fout.write("\\t".join(r1) + "\\n")'
		pSimRead.args.delimit = ["\t", "|"]
		PyPPL(config).start(pSimRead).run()
		procOK(pSimRead, 'simread.txt', self)

if __name__ == '__main__':
	unittest.main()
		