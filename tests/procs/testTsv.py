import helpers, unittest, testly
from os import path
from pyppl import PyPPL
from bioprocs.tsv import pMatrixR, pCbind, pRbind, pCsplit, pRsplit, pSimRead, pTsv
from helpers import getfile, procOK, config

class TestTsv (helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestTsv')
	
	def testMatrix(self):
		name               = 'matrix.txt'
		pMatrixR1           = pMatrixR.copy()
		pMatrixR1.input     = [getfile(name)]
		pMatrixR1.args.code = "mat = log(mat + 1, 2)"
		PyPPL(config).start(pMatrixR1).run()
		procOK(pMatrixR1, name, self)

	def testMatrixNoHead(self):
		name                 = 'matrix-nohead.txt'
		pMatrixR2             = pMatrixR.copy()
		pMatrixR2.input       = [getfile(name)]
		pMatrixR2.args.code   = "mat = log(mat + 1, 2)"
		pMatrixR2.args.inopts.cnames = False
		PyPPL(config).start(pMatrixR2).run()
		self.assertFileEqual(pMatrixR2.channel.get(), path.join(self.outdir, name))

	def testMatrixNoRn(self):
		name                 = 'matrix-norn.txt'
		pMatrixR3             = pMatrixR.copy()
		pMatrixR3.input       = [getfile(name)]
		pMatrixR3.args.code   = "mat = log(mat + 1, 2)"
		pMatrixR3.args.inopts.rnames = False
		PyPPL(config).start(pMatrixR3).run()
		self.assertFileEqual(pMatrixR3.channel.get(), path.join(self.outdir, name))

	def testCbind(self):
		name1 = 'matrix-cbind-1.txt'
		name2 = 'matrix-cbind-2.txt'
		name0 = 'matrix-cbind.txt'
		pCbind1 = pCbind.copy()
		pCbind1.input = [[getfile(name1), getfile(name2)]]
		PyPPL(config).start(pCbind1).run()
		self.assertFileEqual(pCbind1.channel.get(), path.join(self.outdir, name0))

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
		procOK(pCsplit, 'matrix-rbind-1.csplits', self)

	def testRsplit(self):
		name1 = 'matrix-rbind-2.txt'
		pRsplit.input = [getfile(name1)]
		PyPPL(config).start(pRsplit).run()
		procOK(pRsplit, 'matrix-rbind-2.rsplits', self)

	def testRsplitN(self):
		name1 = 'matrix-rsplit.txt'
		pRsplitN = pRsplit.copy()
		pRsplitN.input = [getfile(name1)]
		pRsplitN.args.size = 3
		PyPPL(config).start(pRsplitN).run()
		procOK(pRsplitN, 'matrix-rsplit.rsplits', self)

	def testCsplitN(self):
		name1 = 'matrix-csplit.txt'
		pCsplitN = pCsplit.copy()
		pCsplitN.input = [getfile(name1)]
		pCsplitN.args.size = 3
		PyPPL(config).start(pCsplitN).run()
		procOK(pCsplitN, 'matrix-csplit.csplits', self)
		
	def dataProvider_testSimRead(self):
		infile1 = path.join(self.indir, 'simread1.txt.gz')
		infile2 = path.join(self.indir, 'simread2.txt')
		exptfile = path.join(self.outdir, 'simread.txt')
		args = {
			'inopts': {'skip': [3], 'delimit': ["\t", "|"]},
			'usemeta': 1,
			'match': 'lambda r1, r2: SimRead.compare(r1[0], r2[1])',
			'do': 'lambda r1, r2: writer.write(r1)'
		}
		yield 't1', infile1, infile2, args, exptfile

	def testSimRead(self, tag, infile1, infile2, args, exptfile):
		pSimReadTest = pSimRead.copy(tag = tag)
		pSimReadTest.input = [[infile1, infile2]]
		pSimReadTest.args.update(args)
		PyPPL(config).start(pSimReadTest).run()
		self.assertFileEqual(pSimReadTest.channel.get(), exptfile)
	
	def dataProvider_testTsv(self):
		infile = path.join(self.indir, 'matrix-nohead.txt')
		exptfile = path.join(self.outdir, 'matrix-nohead.tsv')
		args = {
			'inopts': {'cnames': ['NAME', 'COL1', 'COL2', 'COL3']},
			'outopts': {'head': True}
		}
		yield 't1', infile, args, exptfile
	
	def testTsv(self, tag, infile, args, exptfile):
		pTsvTest = pTsv.copy(tag = tag)
		pTsvTest.input = [infile]
		pTsvTest.args.update(args)
		PyPPL(config).start(pTsvTest).run()
		self.assertFileEqual(pTsvTest.channel.get(), exptfile)


if __name__ == '__main__':
	testly.main(failfast = True)
		