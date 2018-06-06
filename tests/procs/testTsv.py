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

	def dataProvider_testCbind(self):
		infiles1 = [path.join(self.indir, infile) for infile in ['matrix-cbind-1.txt', 'matrix-cbind-2.txt']]
		outfile1 = path.join(self.outdir, 'matrix-cbind.txt')
		outfile2 = path.join(self.outdir, 'matrix-cbind2.txt')
		infiles2 = [path.join(self.indir, infile) for infile in ['matrix-cbind-1.txt', 'matrix-cbind-3.txt']]
		infiles2 = [path.join(self.indir, infile) for infile in ['matrix-cbind-1.txt', 'matrix-cbind-3.txt']]
		infiles3 = [path.join(self.indir, infile) for infile in ['matrix-cbind-1.txt', 'matrix-cbind-4.txt']]
		outfile3 = path.join(self.outdir, 'matrix-cbind3.txt')
		outfile4 = path.join(self.outdir, 'matrix-cbind4.txt')
		# default with fill
		yield 't1', infiles1, outfile1, {}
		# no rownames and skip the first line
		yield 't2', infiles1, outfile2, {'fill': True, 'inopts': {'rnames': False, 'skip': 1}}
		# different delimit
		yield 't3', infiles2, outfile2, {'fill': True, 'inopts': {'rnames': False, 'skip': 1, 'delimit': ['\t', ';']}}
		# no fill with rownames (at the same order)
		yield 't4', infiles3, outfile3, {'fill': False, 'inopts': {'rnames': True}}
		# use file name as column name
		yield 't5', infiles1, outfile4, {'fill': True, 'inopts': {'cnames': False, 'skip': 1}}

	def testCbind(self, tag, infiles, outfile, args):
		pCbindTest = pCbind.copy(tag = tag)
		if 'inopts' in args:
			pCbindTest.args.inopts.update(args['inopts'])
			del args['inopts']
		if 'params' in args:
			pCbindTest.args.params.update(args['params'])
			del args['params']
		pCbindTest.args.update(args)
		pCbindTest.input = [infiles]
		PyPPL(config).start(pCbindTest).run()
		self.assertFileEqual(pCbindTest.channel.get(), outfile)

	def dataProvider_testRbind(self):
		infiles1 = [path.join(self.indir, infile) for infile in ['matrix-rbind-1.txt', 'matrix-rbind-2.txt']]
		outfile1 = path.join(self.outdir, 'matrix-rbind.txt')
		outfile2 = path.join(self.outdir, 'matrix-rbind2.txt')
		infiles2 = [path.join(self.indir, infile) for infile in ['matrix-rbind-1.txt', 'matrix-rbind-3.txt']]
		infiles2 = [path.join(self.indir, infile) for infile in ['matrix-rbind-1.txt', 'matrix-rbind-3.txt']]
		infiles3 = [path.join(self.indir, infile) for infile in ['matrix-rbind-1.txt', 'matrix-rbind-4.txt']]
		outfile3 = path.join(self.outdir, 'matrix-rbind3.txt')
		outfile4 = path.join(self.outdir, 'matrix-rbind4.txt')
		# default with fill
		yield 't1', infiles1, outfile1, {}
		# no colnames and skip the first line
		yield 't2', infiles1, outfile2, {'fill': True, 'inopts': {'cnames': False, 'skip': 1}}
		# different delimit
		yield 't3', infiles2, outfile2, {'fill': True, 'inopts': {'cnames': False, 'skip': 1, 'delimit': ['\t', ';']}}
		# no fill with colnames (at the same order)
		yield 't4', infiles3, outfile3, {'fill': False, 'inopts': {'cnames': True}}
		yield 't5', infiles1, outfile4, {'fill': True, 'inopts': {'rnames': False}}

	def testRbind(self, tag, infiles, outfile, args):
		pRbindTest = pRbind.copy(tag = tag)
		if 'inopts' in args:
			pRbindTest.args.inopts.update(args['inopts'])
			del args['inopts']
		if 'params' in args:
			pRbindTest.args.params.update(args['params'])
			del args['params']
		pRbindTest.args.update(args)
		pRbindTest.input = [infiles]
		PyPPL(config).start(pRbindTest).run()
		self.assertFileEqual(pRbindTest.channel.get(), outfile)

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
		